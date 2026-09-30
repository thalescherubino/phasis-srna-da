"""DA from a completed, Phasis-backed CLI count run; never alter source results."""

from __future__ import annotations

import csv
import json
import sys
from importlib.metadata import version
from pathlib import Path

from .errors import ValidationError
from .io import read_targets
from .provenance import base_run_metadata, utc_timestamp, write_input_manifest, write_json


REQUIRED_PROJECTIONS = {
    "global_atomic", "phasi_product", "phasi_locus", "phasi_locus_tag", "phasi_tag_seq"
}
INHERITED_FIELDS = (
    "assembly_id", "annotation_id", "library_format", "library_min_length",
    "library_max_length", "ambiguous_N_policy", "reference_orientation", "mapping_orientation",
)


def _table(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise ValidationError(f"Missing count-run evidence: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _integer(value: str) -> int:
    try:
        number = int(value)
    except (ValueError, TypeError) as exc:
        raise ValidationError("Count-run counts and depths must be nonnegative integers.") from exc
    if number < 0:
        raise ValidationError("Count-run counts and depths must be nonnegative integers.")
    return number


def check_count_run(source: Path):
    """Check lineage, complete sample coverage, and raw matrix conservation before fitting."""

    metadata_path = source / "provenance/run_metadata.json"
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ValidationError(f"Cannot read count-run metadata: {metadata_path}") from exc
    if not isinstance(metadata, dict) or metadata.get("tool") != "phasis-srna-da" or metadata.get("status") != "complete":
        raise ValidationError("--count-run must be a completed phasis-srna-da count run.")
    if any(metadata.get(key) is None for key in INHERITED_FIELDS):
        raise ValidationError("Count run is missing explicit reference/length/filtering provenance.")
    inputs_path = source / "provenance/input_manifest.tsv"
    phasis = [row for row in _table(inputs_path) if row.get("role") == "phasis_product_catalog"]
    if not phasis or any(not r.get("path") or len(r.get("sha256", "")) != 64 for r in phasis):
        raise ValidationError("Count run must retain checksummed mandatory Phasis product input provenance.")
    targets = source / "provenance/resolved_targets.tsv"
    if not targets.is_file():
        raise ValidationError(f"Missing count-run targets: {targets}")
    samples = tuple(read_targets(targets))
    sample_ids = {sample.sample_id for sample in samples}
    summary_path = source / "qc/library_summary.tsv"
    qc_rows = _table(summary_path)
    qc = {r["sample_id"]: r for r in qc_rows}
    if set(qc) != sample_ids or len(qc_rows) != len(qc):
        raise ValidationError("Count-run QC does not contain exactly the declared samples.")
    for sample in samples:
        row = qc[sample.sample_id]
        retained = _integer(row["retained_abundance"])
        if retained < 1 or _integer(row["normalization_depth"]) != retained:
            raise ValidationError(f"{sample.sample_id}: normalization depth must equal positive retained abundance.")
        if row["condition"] != sample.condition or row["biological_replicate"] != sample.biological_replicate:
            raise ValidationError(f"{sample.sample_id}: source targets and QC design disagree.")
        if sum(_integer(row[k]) for k in (
            "primary_unique_assigned_abundance", "primary_equivalence_assigned_abundance",
            "primary_unassigned_abundance",
        )) != retained:
            raise ValidationError(f"{sample.sample_id}: QC abundance conservation failed.")

    matrix_paths = sorted((source / "counts").glob("*.raw.tsv"))
    projections = {p.name.removesuffix(".raw.tsv") for p in matrix_paths}
    if not REQUIRED_PROJECTIONS.issubset(projections):
        raise ValidationError("Count run is missing mandatory global/Phasis count projections.")
    conservation_path = source / "mapping/conservation_by_sample.tsv"
    rows = _table(conservation_path)
    cons = {(r["projection"], r["sample_id"]): r for r in rows}
    if set(cons) != {(p, s) for p in projections for s in sample_ids} or len(rows) != len(cons):
        raise ValidationError("Count-run conservation evidence has incomplete or duplicate sample/projection rows.")
    for path in matrix_paths:
        projection = path.name.removesuffix(".raw.tsv")
        totals = dict.fromkeys(sample_ids, 0)
        seen = set()
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not (sample_ids | {"unit_id"}).issubset(reader.fieldnames or []):
                raise ValidationError(f"{path}: count matrix is missing unit_id or target samples.")
            for row in reader:
                unit = row["unit_id"]
                if not unit or unit in seen:
                    raise ValidationError(f"{path}: blank or duplicate count unit.")
                seen.add(unit)
                for sample in sample_ids:
                    totals[sample] += _integer(row[sample])
        for sample, total in totals.items():
            evidence = cons[projection, sample]
            retained = _integer(qc[sample]["retained_abundance"])
            assigned = sum(_integer(evidence[k]) for k in ("unique_assigned_abundance", "equivalence_assigned_abundance"))
            if not (
                total == assigned == _integer(evidence["assigned_abundance"])
                and assigned + _integer(evidence["unassigned_abundance"]) == retained
                and _integer(evidence["input_abundance"]) == retained
                and evidence["conserved"] == "true"
            ):
                raise ValidationError(f"{projection}/{sample}: count matrix abundance conservation failed.")
            if projection == "global_atomic" and any(
                _integer(evidence[k]) != _integer(qc[sample]["primary_" + k])
                for k in ("unique_assigned_abundance", "equivalence_assigned_abundance", "unassigned_abundance")
            ):
                raise ValidationError(f"{sample}: global matrix and library QC disagree.")
    paths = [metadata_path, inputs_path, targets, summary_path, conservation_path, *matrix_paths]
    return metadata, samples, paths


def command_da(args) -> int:
    from .cli import _validate_declared_contrasts
    from .da import _optional_dependencies, run_differential_abundance

    if not 0 < args.fdr < 1:
        raise ValidationError("--fdr must be between 0 and 1.")
    if args.outdir == args.count_run or args.outdir.is_relative_to(args.count_run):
        raise ValidationError("DA output must be separate from the read-only count run.")
    if args.outdir.exists():
        raise FileExistsError(f"DA output already exists: {args.outdir}; choose a new directory.")
    try:
        source_metadata, samples, evidence_paths = check_count_run(args.count_run)
    except KeyError as exc:
        raise ValidationError(f"Count-run evidence is missing required field: {exc}") from exc
    _validate_declared_contrasts(args.contrast, {s.condition for s in samples})
    _optional_dependencies()
    metadata = base_run_metadata(command_line=getattr(args, "_command_line", list(sys.argv)))
    metadata.update({key: source_metadata[key] for key in INHERITED_FIELDS})
    metadata.update({
        "analysis_scope": "DA from completed mandatory-Phasis-backed count run; no remapping or discovery",
        "source_count_run": str(args.count_run),
        "da_normalization_policy": args.normalization,
        "da_contrasts": [{"numerator": n, "denominator": d} for n, d in args.contrast],
        "min_count": args.min_count, "prefilter_policy": args.prefilter_policy,
        "fdr": args.fdr, "threads": args.threads,
        "software_versions": {name: version(name) for name in ("pydeseq2", "numpy", "pandas")},
        "status": "running",
    })
    metadata_path = args.outdir / "provenance/run_metadata.json"
    write_json(metadata_path, metadata)
    try:
        write_input_manifest(args.outdir, [
            ("source_count_evidence", path, "Read-only source; preserve mandatory Phasis lineage and retained-depth normalization")
            for path in evidence_paths
        ])
        run_differential_abundance(
            count_directory=args.count_run / "counts", samples=samples, outdir=args.outdir / "da",
            min_count=args.min_count, alpha=args.fdr, threads=args.threads,
            library_summary_path=args.count_run / "qc/library_summary.tsv", contrasts=tuple(args.contrast),
            prefilter_policy=args.prefilter_policy,
        )
    except Exception as exc:
        metadata.update(status="failed", failed_utc=utc_timestamp(), error=str(exc))
        write_json(metadata_path, metadata)
        raise
    metadata.update(status="complete", completed_utc=utc_timestamp())
    write_json(metadata_path, metadata)
    print(f"Completed declared DA comparisons: {args.outdir}")
    return 0
