"""Local command-line entry point. No network or genomic uploads."""
import argparse
import json
from pathlib import Path
import sys
import tempfile
import time

from . import __version__
from .models import PhaseSpaceError


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="phasespace", description="Explore possible coding consequences of called variants across unresolved phase blocks.")
    parser.add_argument("--version", action="version", version=f"PhaseSpace {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)
    report = commands.add_parser("report", help="Analyze a local VCF/gVCF locus and write an interactive HTML report")
    for name, help_text in [("vcf", "VCF, bgzip VCF/gVCF, or BCF"), ("fasta", "Matching reference FASTA with .fai index"), ("annotation", "GTF or GFF3 annotation"), ("transcript", "Exact transcript ID")]:
        report.add_argument(f"--{name}", required=True, help=help_text)
    report.add_argument("--sample", help="Sample name (required for multisample files)")
    report.add_argument("--min-gq", type=float, default=20, help="Minimum genotype quality; missing values fail (default: 20)")
    report.add_argument("--min-dp", type=int, default=5, help="Minimum depth, using MIN_DP for blocks when present (default: 5)")
    report.add_argument("--phase-from", metavar="PHASED_VCF", help="Phased VCF/BCF (e.g. from WhatsHap or longphase) that supplies GT order and PS for heterozygous calls; the main input still supplies evidence")
    report.add_argument("--max-blocks", type=int, default=8, help="Hard cap on independent heterozygous phase blocks, 0–12 (default: 8)")
    demo = commands.add_parser("demo", help="Generate a report from a synthetic two-variant locus")
    demo.add_argument("--write-inputs", type=Path, help="Keep synthetic FASTA, GTF, VCF and indexed gVCF in this new directory")
    for command in (report, demo):
        command.add_argument("--out", required=True, type=Path, help="Output standalone HTML path")
        command.add_argument("--json", type=Path, help="Also write complete machine-readable results")
        command.add_argument("--force", action="store_true", help="Replace existing output files")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    try:
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
            sources = [p for p in (args.vcf, args.fasta, args.annotation, args.phase_from) if p is not None]
            inputs = {Path(p).resolve() for p in sources}
            inputs.update(Path(str(p) + suffix).resolve() for p in sources for suffix in (".fai", ".gzi", ".tbi", ".csi"))
            if any(path.resolve() in inputs for path in outputs):
                raise PhaseSpaceError("Output paths must not overwrite input files.")
            locus = load_locus(args.vcf, args.fasta, args.annotation, args.transcript, sample=args.sample, min_gq=args.min_gq, min_dp=args.min_dp, phase_from=args.phase_from)
            result = analyze(locus, max_blocks=args.max_blocks)
            result["synthetic"] = False
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
