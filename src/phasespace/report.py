"""Write a portable, self-contained report without remote assets or requests."""
from __future__ import annotations

import json
from importlib.resources import files
from pathlib import Path


def _script_json(value: dict) -> str:
    """JSON embedded in HTML must not contain a script-closing delimiter."""
    return (
        json.dumps(value, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
        .replace("&", "\\u0026")
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("\u2028", "\\u2028")
        .replace("\u2029", "\\u2029")
    )


def _compact(result: dict) -> dict:
    """Drop per-copy sequences that the browser can rebuild from called edits.

    Two full sequences per arrangement grow as 2**blocks x CDS length (tens of
    megabytes at the block cap). Each copy keeps its variant IDs and chosen
    alleles, so report.js reconstructs sequences exactly on demand.
    """
    scenarios = result.get("scenarios", [])
    copies = [copy for scenario in scenarios for copy in scenario.get("haplotypes", [])]
    if not copies or any("alleles" not in copy for copy in copies):
        return result
    compact = dict(result)
    compact["scenarios"] = [
        dict(scenario, haplotypes=[
            {key: value for key, value in copy.items() if key not in ("cds", "protein")}
            for copy in scenario["haplotypes"]
        ])
        for scenario in scenarios
    ]
    compact["report_encoding"] = "edits"
    return compact


def write_report(result: dict, output_path: str | Path) -> None:
    """Render an analysis dictionary to a standalone interactive HTML file.

    Input strings are only embedded as inert JSON. The browser renderer inserts
    every data-derived label as text, including variant IDs and sample names.
    """
    assets = files("phasespace").joinpath("assets")
    template = assets.joinpath("report.html").read_text(encoding="utf-8")
    css = assets.joinpath("report.css").read_text(encoding="utf-8")
    javascript = assets.joinpath("report.js").read_text(encoding="utf-8")
    # Substitute data last: strings inside input metadata cannot become template
    # substitutions, even if they happen to contain one of these markers.
    html = template.replace("/* PHASESPACE_CSS */", css)
    html = html.replace("/* PHASESPACE_JS */", javascript)
    html = html.replace("PHASESPACE_DATA", _script_json(_compact(result)))
    Path(output_path).write_text(html, encoding="utf-8")
