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
    html = html.replace("PHASESPACE_DATA", _script_json(result))
    Path(output_path).write_text(html, encoding="utf-8")
