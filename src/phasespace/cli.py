"""Local command-line entry point. No network or genomic uploads."""
import argparse
import json
from pathlib import Path
import sys
import tempfile
import time

from . import __version__
from .models import PhaseSpaceError


def _inputs(command: argparse.ArgumentParser) -> None:
    command.add_argument("--vcf", required=True, help="VCF, bgzip VCF/gVCF, or BCF")
    command.add_argument("--fasta", required=True, help="Matching reference FASTA with .fai index")
    command.add_argument("--annotation", required=True, help="GTF or GFF3 annotation")
    command.add_argument("--sample", help="Sample name (required for multisample files)")
    command.add_argument("--min-gq", type=float, default=20, help="Minimum genotype quality; missing values fail (default: 20)")
    command.add_argument("--min-dp", type=int, default=5, help="Minimum depth, using MIN_DP for blocks when present (default: 5)")
    command.add_argument("--phase-from", metavar="PHASED_VCF", help="Phased VCF/BCF (e.g. from WhatsHap or longphase) that supplies GT order and PS for heterozygous calls; the main input still supplies evidence")
    command.add_argument("--max-blocks", type=int, default=8, help="Hard cap on independent heterozygous phase blocks, 0–12 (default: 8)")
    command.add_argument("--reads", metavar="BAM", action="append", default=[], help="Indexed BAM/CRAM; counts reads that show which alleles travel together (repeatable)")
    command.add_argument("--info", metavar="FIELDS", help="Comma-separated INFO fields to show for each variant (e.g. CLNSIG,AF)")
    command.add_argument("--info-from", metavar="ANNOTATED_VCF", help="Indexed VCF (e.g. ClinVar or gnomAD sites) to read --info fields from, matched by position, REF and ALT; default: the main VCF")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="phasespace", description="Explore possible coding consequences of called variants across unresolved phase blocks.")
    parser.add_argument("--version", action="version", version=f"PhaseSpace {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)
    report = commands.add_parser("report", help="Analyze one transcript and write an interactive HTML report")
    _inputs(report)
    which = report.add_mutually_exclusive_group(required=True)
    which.add_argument("--transcript", help="Transcript ID as written in the annotation, with or without version")
    which.add_argument("--gene", help="Gene symbol or ID; uses its MANE Select, else Ensembl canonical, else RefSeq Select transcript")
    demo = commands.add_parser("demo", help="Generate a report from a synthetic two-variant locus")
    demo.add_argument("--write-inputs", type=Path, help="Keep synthetic FASTA, GTF, VCF and indexed gVCF in this new directory")
    for command in (report, demo):
        command.add_argument("--out", required=True, type=Path, help="Output standalone HTML path")
        command.add_argument("--json", type=Path, help="Also write complete machine-readable results")
        command.add_argument("--force", action="store_true", help="Replace existing output files")
    scan = commands.add_parser("scan", help="Scan every gene (or a gene list) and rank where phase matters")
    _inputs(scan)
    scan.add_argument("--genes", help="Comma-separated gene symbols/IDs, or @FILE with one per line (default: all protein-coding genes)")
    scan.add_argument("--out-dir", required=True, type=Path, help="Directory for index.html, scan.tsv, scan.json and reports/")
    scan.add_argument("--reports", choices=("candidates", "all", "none"), default="candidates",
                      help="Write reports for candidate genes (phase decides, or no reference protein with a disruptive change), all genes with variants, or none (default: candidates)")
    scan.add_argument("--threads", type=int, default=None, help="Worker processes (default: CPUs − 1, at most 8)")
    scan.add_argument("--force", action="store_true", help="Write into an existing non-empty output directory")
    return parser


def _gene_list(value: str | None) -> set[str] | None:
    if not value:
        return None
    if value.startswith("@"):
        text = Path(value[1:]).read_text(encoding="utf-8")
        names = [line.split("#", 1)[0].strip() for line in text.splitlines()]
    else:
        names = [name.strip() for name in value.split(",")]
    genes = {name for name in names if name}
    if not genes:
        raise PhaseSpaceError("--genes is empty.")
    return genes


def _fields(value: str | None) -> tuple[str, ...]:
    return tuple(name.strip() for name in (value or "").split(",") if name.strip())


def _check_not_inputs(args: argparse.Namespace, outputs: list[Path]) -> None:
    sources = [p for p in (args.vcf, args.fasta, args.annotation, args.phase_from, args.info_from, *args.reads) if p is not None]
    inputs = {Path(p).resolve() for p in sources}
    inputs.update(Path(str(p) + suffix).resolve() for p in sources for suffix in (".fai", ".gzi", ".tbi", ".csi"))
    if any(path.resolve() in inputs for path in outputs):
        raise PhaseSpaceError("Output paths must not overwrite input files.")


def _scan(args: argparse.Namespace) -> int:
    from .analysis import check_max_blocks
    from .inputs import check_thresholds
    from .scan import default_threads, run_scan
    check_max_blocks(args.max_blocks)
    check_thresholds(args.min_gq, args.min_dp)
    if args.out_dir.exists() and any(args.out_dir.iterdir()) and not args.force:
        raise PhaseSpaceError(f"Output directory is not empty: {args.out_dir}. Use --force to write into it.")
    _check_not_inputs(args, [args.out_dir / "index.html", args.out_dir / "scan.tsv", args.out_dir / "scan.json"])
    threads = args.threads if args.threads is not None else default_threads()
    if threads < 1:
        raise PhaseSpaceError("--threads must be at least 1.")
    result = run_scan(args.vcf, args.fasta, args.annotation, args.out_dir, sample=args.sample, genes=_gene_list(args.genes),
                      phase_from=args.phase_from, info_fields=_fields(args.info), info_from=args.info_from, reads=args.reads,
                      min_gq=args.min_gq, min_dp=args.min_dp, max_blocks=args.max_blocks,
                      reports=args.reports, threads=threads)
    counts = result["counts"]
    print(f"PhaseSpace scan: {counts['transcripts_scanned']:,} transcripts, {counts['with_coding_variants']:,} with coding variants")
    print(f"  phase decides: {counts['phase_decides']:,}   no reference protein (disruptive): {counts['no_reference_protein_disruptive']:,}   not analysed: {counts['not_analysed']:,}")
    print(f"Index: {(args.out_dir / 'index.html').resolve()}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "scan":
            return _scan(args)
        from .analysis import analyze, check_max_blocks
        from .inputs import check_thresholds, load_locus
        from .report import write_report
        outputs = [p for p in (args.out, args.json) if p is not None]
        if len({p.resolve() for p in outputs}) != len(outputs):
            raise PhaseSpaceError("HTML and JSON output paths must be different.")
        for path in outputs:
            if path.exists() and not args.force:
                raise PhaseSpaceError(f"Output already exists: {path}. Use --force to replace it.")
        started = time.perf_counter()
        if args.command == "demo":
            from .demo import write_demo
            with tempfile.TemporaryDirectory(prefix="phasespace-demo-") as temp:
                files = write_demo(args.write_inputs or temp)
                locus = load_locus(str(files["vcf"]), str(files["fasta"]), str(files["annotation"]), "DEMO1.1", sample="SYNTHETIC")
                locus.source = "Synthetic PhaseSpace demonstration"
                result = analyze(locus)
            result["synthetic"] = True
        else:
            # Validate before reading any input files.
            check_max_blocks(args.max_blocks)
            check_thresholds(args.min_gq, args.min_dp)
            _check_not_inputs(args, outputs)
            transcript, selection = args.transcript, None
            if args.gene:
                from .annotation import select_transcript
                transcript, selection = select_transcript(Path(args.annotation), args.gene)
            locus = load_locus(args.vcf, args.fasta, args.annotation, transcript, sample=args.sample, min_gq=args.min_gq, min_dp=args.min_dp,
                               phase_from=args.phase_from, info_fields=_fields(args.info), info_from=args.info_from)
            if selection:
                locus.warnings.insert(0, selection)
            result = analyze(locus, max_blocks=args.max_blocks)
            result["synthetic"] = False
            if args.reads:
                from .reads import add_read_evidence
                add_read_evidence(result, args.reads, args.fasta)
        result["version"] = __version__
        result["analysis_seconds"] = round(time.perf_counter() - started, 6)
        for path in outputs:
            path.parent.mkdir(parents=True, exist_ok=True)
        write_report(result, args.out)
        if args.json:
            args.json.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
        print(f"PhaseSpace: {result['scenario_count']} arrangements, {result['outcome_count']} distinct coding outcomes")
        print(f"Report: {args.out.resolve()}")
        print(f"Analysis: {result['analysis_seconds']:.3f} s (local inputs)")
        for warning in result.get("warnings", []):
            print(f"Note: {warning}", file=sys.stderr)
        return 0
    except (PhaseSpaceError, OSError) as exc:
        print(f"phasespace: {exc}", file=sys.stderr)
        return 2
