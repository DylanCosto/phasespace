"""Scan every gene (or a gene list) and rank where phase matters.

One representative transcript per protein-coding gene (MANE Select, then
Ensembl canonical, then RefSeq Select) is analysed with the same model as a
single report. Output: ``scan.tsv`` (genes with called coding variants or
problems), ``index.html`` (sortable, linking to reports) and ``reports/``.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
import csv
from importlib.resources import files
import json
import os
from pathlib import Path
import re
import sys
import time

from . import __version__
from .analysis import analyze, check_max_blocks
from .annotation import TranscriptModel, catalogue
from .inputs import LocusReader, _resolve_contig, check_thresholds
from .models import PhaseSpaceError
from .report import _script_json, write_report

TIERS = {1: "Phase decides", 2: "No reference protein, disruptive", 3: "No reference protein", 4: "Reference protein kept"}
COLUMNS = ["tier", "gene", "transcript", "selection", "contig", "verdict", "variants", "heterozygous", "phased",
           "blocks", "arrangements", "protein_outcomes", "changes", "disruptive", "splice_region", "info", "unknown_cds_bases",
           "cds_bases", "status", "report"]

_worker: dict = {}


def _start_worker(settings: dict) -> None:
    _worker["settings"] = settings
    _worker["reader"] = LocusReader(settings["vcf"], settings["fasta"], settings["sample"], settings["min_gq"],
                                    settings["min_dp"], settings["phase_from"], settings["info_fields"], settings["info_from"])


def _safe(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)[:80]


def _scan_one(model: TranscriptModel) -> dict:
    settings = _worker["settings"]
    row = {"gene": model.gene, "transcript": model.id, "selection": model.selection, "contig": model.contig,
           "cds_bases": sum(end - start for start, end in model.intervals), "report": "", "status": "ok"}
    try:
        locus = _worker["reader"].load(model)
    except PhaseSpaceError as exc:
        return {**row, "tier": 5, "status": f"not analysed: {exc}"}
    variants = locus.variants
    row.update(contig=locus.transcript.contig, variants=len(variants),
               heterozygous=sum(v.gt[0] != v.gt[1] for v in variants),
               phased=sum(v.phased and v.gt[0] != v.gt[1] for v in variants),
               unknown_cds_bases=sum(e.end - e.start for e in locus.evidence if e.status == "unknown"),
               splice_region=len(locus.splice_region),
               info="; ".join(f"{v.id.split(':', 1)[1]} {k}={value}" for v in variants for k, value in v.info.items()))
    if not variants:
        return {**row, "tier": 4, "verdict": "good", "status": "no called coding variants"}
    try:
        result = analyze(locus, max_blocks=settings["max_blocks"])
    except PhaseSpaceError as exc:
        return {**row, "tier": 5, "status": f"not analysed: {exc}"}
    summary = result["summary"]
    row.update(tier=summary["tier"], verdict=summary["verdict"], blocks=len(result["blocks"]),
               arrangements=result["scenario_count"], protein_outcomes=result["protein_outcome_count"],
               changes=",".join(summary["changes"]), disruptive=",".join(summary["disruptive"]))
    wanted = settings["reports"] == "all" or (settings["reports"] == "candidates" and summary["tier"] <= 2)
    if wanted:
        result.update(synthetic=False, version=__version__, analysis_seconds=0.0)
        name = f"reports/{_safe(model.gene)}_{_safe(model.id)}.html"
        write_report(result, Path(settings["out_dir"]) / name)
        row["report"] = name
    return row


def run_scan(vcf: str, fasta: str, annotation: str, out_dir: str | Path, *, sample: str | None = None,
             genes: set[str] | None = None, phase_from: str | None = None, info_fields: tuple[str, ...] = (),
             info_from: str | None = None, min_gq: float = 20, min_dp: int = 5,
             max_blocks: int = 8, reports: str = "candidates", threads: int = 1, log=None) -> dict:
    check_thresholds(min_gq, min_dp)
    check_max_blocks(max_blocks)
    if reports not in ("candidates", "all", "none"):
        raise PhaseSpaceError("--reports must be candidates, all or none.")
    out_dir = Path(out_dir)
    (out_dir / "reports").mkdir(parents=True, exist_ok=True)
    log = log or (lambda message: print(message, file=sys.stderr))
    started = time.perf_counter()
    log("Reading annotation …")
    models, skipped = catalogue(Path(annotation), genes)
    all_models = models
    if not models:
        raise PhaseSpaceError("No supported coding transcripts were found to scan.")
    settings = {"vcf": str(vcf), "fasta": str(fasta), "sample": sample, "min_gq": min_gq, "min_dp": min_dp,
                "phase_from": phase_from, "info_fields": tuple(info_fields), "info_from": info_from, "max_blocks": max_blocks, "reports": reports, "out_dir": str(out_dir)}
    # Fail fast on unreadable inputs or an ambiguous sample before starting workers.
    with LocusReader(vcf, fasta, sample, min_gq, min_dp, phase_from, info_fields, info_from) as probe:
        settings["sample"] = sample = probe.sample
        source = probe.source
        references = set(probe.fasta.references)
        declared = set(probe.vcf.header.contigs)
    # Inputs limited to some contigs (one chromosome, an exome region set) are
    # normal; genes elsewhere are counted, not reported as failures.
    present = [m for m in models if _resolve_contig(m.contig, references) and (not declared or _resolve_contig(m.contig, declared))]
    absent_contigs = sorted({m.contig for m in models} - {m.contig for m in present})
    models = present
    skipped = [entry for entry in skipped if not entry[2] or (_resolve_contig(entry[2], references) and (not declared or _resolve_contig(entry[2], declared)))]
    if not models:
        raise PhaseSpaceError("None of the annotated contigs are present in both the reference FASTA and the VCF.")
    log(f"Scanning {len(models):,} transcripts with {threads} process{'es' if threads != 1 else ''} …")
    rows: list[dict] = []
    if threads <= 1:
        _start_worker(settings)
        try:
            for index, model in enumerate(models, 1):
                rows.append(_scan_one(model))
                if index % 2000 == 0:
                    log(f"  {index:,} / {len(models):,}")
        finally:
            _worker["reader"].close()
    else:
        with ProcessPoolExecutor(threads, initializer=_start_worker, initargs=(settings,)) as pool:
            for index, row in enumerate(pool.map(_scan_one, models, chunksize=16), 1):
                rows.append(row)
                if index % 2000 == 0:
                    log(f"  {index:,} / {len(models):,}")
    rows.sort(key=lambda row: (row["tier"], -len(row.get("disruptive", "") or ""), row["gene"]))
    shown = [row for row in rows if row.get("variants") or row.get("splice_region") or row["tier"] == 5]
    with open(out_dir / "scan.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, COLUMNS, delimiter="\t", extrasaction="ignore", restval="")
        writer.writeheader()
        writer.writerows(shown)
    counts = {
        "transcripts_scanned": len(rows),
        "with_coding_variants": sum(1 for row in rows if row.get("variants")),
        "phase_decides": sum(1 for row in rows if row["tier"] == 1),
        "no_reference_protein_disruptive": sum(1 for row in rows if row["tier"] == 2),
        "no_reference_protein": sum(1 for row in rows if row["tier"] == 3),
        "not_analysed": sum(1 for row in rows if row["tier"] == 5),
        "genes_skipped": len(skipped),
        "genes_on_absent_contigs": len(all_models) - len(models),
        "reports": sum(1 for row in rows if row["report"]),
    }
    summary = {"sample": sample, "source": source, "annotation": Path(annotation).name, "version": __version__,
               "counts": counts, "tiers": TIERS, "seconds": round(time.perf_counter() - started, 1),
               "skipped": [{"gene": gene, "reason": reason} for gene, reason, _ in skipped][:500],
               "absent_contigs": absent_contigs}
    _write_index(out_dir / "index.html", summary, shown)
    (out_dir / "scan.json").write_text(json.dumps({**summary, "rows": shown}, indent=1) + "\n", encoding="utf-8")
    return {**summary, "rows": shown}


def _write_index(path: Path, summary: dict, rows: list[dict]) -> None:
    assets = files("phasespace").joinpath("assets")
    html = assets.joinpath("scan.html").read_text(encoding="utf-8")
    html = html.replace("/* PHASESPACE_CSS */", assets.joinpath("report.css").read_text(encoding="utf-8") + assets.joinpath("scan.css").read_text(encoding="utf-8"))
    html = html.replace("/* PHASESPACE_JS */", assets.joinpath("scan.js").read_text(encoding="utf-8"))
    html = html.replace("PHASESPACE_DATA", _script_json({**summary, "rows": rows}))
    path.write_text(html, encoding="utf-8")


def default_threads() -> int:
    return max(1, min(8, (os.cpu_count() or 2) - 1))
