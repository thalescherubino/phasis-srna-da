"""Reproducibility manifests for a downstream run that never modifies Phasis."""

from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from . import __version__
from .io import sha256_file, write_tsv


def utc_timestamp() -> str:
    """Return an unambiguous UTC timestamp for a run manifest."""

    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def prepare_output_directory(outdir: Path) -> None:
    """Create an empty result directory without overwriting prior analysis."""

    if outdir.exists():
        if any(outdir.iterdir()):
            raise FileExistsError(
                f"Output directory already exists and is not empty: {outdir}. "
                "Choose a new output directory; this tool never overwrites a prior run."
            )
    else:
        outdir.mkdir(parents=True)
    for name in ("catalogs", "counts", "mapping", "provenance", "qc", "work"):
        (outdir / name).mkdir(exist_ok=True)


def write_json(path: Path, data: Any) -> None:
    """Write sorted, human-readable JSON for a reproducibility record."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8", newline="\n") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
        handle.write("\n")


def base_run_metadata(*, command_line: list[str]) -> dict[str, Any]:
    """Return runtime facts that apply to every output of a CLI run."""

    return {
        "tool": "phasis-srna-da",
        "tool_version": __version__,
        "started_utc": utc_timestamp(),
        "command_line": command_line,
        "python": sys.version,
        "platform": platform.platform(),
        "analysis_scope": (
            "downstream exact re-counting of a user-supplied fixed catalog; "
            "no de novo small-RNA locus discovery"
        ),
    }


def input_manifest_rows(entries: list[tuple[str, Path, str]]) -> list[dict[str, str]]:
    """Checksum every declared input and label its role in a manifest table."""

    rows: list[dict[str, str]] = []
    for role, path, note in entries:
        rows.append(
            {
                "role": role,
                "path": str(path),
                "sha256": sha256_file(path),
                "bytes": str(path.stat().st_size),
                "note": note,
            }
        )
    return rows


def write_input_manifest(outdir: Path, entries: list[tuple[str, Path, str]]) -> None:
    """Write checksummed source paths before a result is considered reusable."""

    write_tsv(
        outdir / "provenance" / "input_manifest.tsv",
        ["role", "path", "sha256", "bytes", "note"],
        input_manifest_rows(entries),
    )
