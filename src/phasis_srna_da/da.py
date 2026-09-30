"""Optional Python-only negative-binomial DA using a pinned PyDESeq2 API."""

from __future__ import annotations

import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from .errors import ValidationError
from .models import Sample


MATRIX_METADATA_COLUMNS = {
    "unit_id",
    "unit_label",
    "assignment_type",
    "projection",
    "srna_class",
    "feature_id",
    "source_kind",
    "phase",
    "identifier",
    "tag_seq",
    "strand",
    "observed_pos",
    "parent_level",
    "parent_id",
}
MIN_FEATURES_FOR_DA = 10
REPORTING_DEPTH = 30_000_000
PREFILTER_ANY_SAMPLES = "any-samples"
PREFILTER_ALL_REPLICATES_ONE_CONDITION = "all-replicates-in-one-condition"
PREFILTER_POLICIES = (
    PREFILTER_ALL_REPLICATES_ONE_CONDITION,
    PREFILTER_ANY_SAMPLES,
)


def _optional_dependencies():
    try:
        import numpy as np
        import pandas as pd
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError as exc:
        raise ValidationError(
            "Optional DA dependencies are unavailable. Install a separate analysis environment with "
            "`python -m pip install -e '.[da]'`; do not alter the Phasis environment."
        ) from exc
    return np, pd, DeseqDataSet, DeseqStats


def _combine_technical_replicates(samples: Iterable[Sample], pd):
    """Return analysis-unit metadata, combining explicitly labelled technical runs.

    A repeated biological replicate without an explicit technical-replicate label
    is rejected: treating it as independent would pseudo-replicate the model.
    """

    groups: dict[tuple[str, str], list[Sample]] = defaultdict(list)
    for sample in samples:
        groups[(sample.condition, sample.biological_replicate)].append(sample)
    mapping: dict[str, str] = {}
    metadata_rows: list[dict[str, str]] = []
    for (condition, biological_replicate), members in sorted(groups.items()):
        if len(members) > 1 and any(not member.technical_replicate for member in members):
            raise ValidationError(
                "Multiple target rows share condition/biological_replicate "
                f"({condition!r}, {biological_replicate!r}) but are not all marked technical_replicate."
            )
        batch_values = {member.batch for member in members if member.batch}
        if len(batch_values) > 1:
            raise ValidationError(
                f"Technical replicates for {condition!r}/{biological_replicate!r} have conflicting batches."
            )
        analysis_id = f"{condition}__bio_rep_{biological_replicate}"
        for member in members:
            mapping[member.sample_id] = analysis_id
        metadata_rows.append(
            {
                "analysis_id": analysis_id,
                "condition": condition,
                "biological_replicate": biological_replicate,
                "batch": next(iter(batch_values), ""),
                "n_technical_libraries": str(len(members)),
            }
        )
    metadata = pd.DataFrame(metadata_rows).set_index("analysis_id")
    return mapping, metadata


def _choose_design(metadata, np, pd) -> tuple[str, str | None]:
    """Use batch only when it is complete and non-confounded with condition."""

    if not metadata["batch"].eq("").all():
        if metadata["batch"].eq("").any():
            return "~condition", "batch ignored because it is missing for at least one analysis sample"
        design_terms = pd.get_dummies(metadata[["batch", "condition"]], drop_first=True, dtype=float)
        design_matrix = np.column_stack([np.ones(len(metadata)), design_terms.to_numpy()])
        if np.linalg.matrix_rank(design_matrix) == design_matrix.shape[1]:
            return "~batch + condition", None
        return "~condition", "batch ignored because it is confounded or rank-deficient with condition"
    return "~condition", None


def _filter_features(
    counts,
    conditions,
    min_count: int,
    policy: str = PREFILTER_ANY_SAMPLES,
):
    """Pre-filter low counts before model fitting, retaining the reason for every row."""

    group_sizes = conditions.value_counts()
    smallest_group = int(group_sizes.min())
    if policy == PREFILTER_ANY_SAMPLES:
        keep = (counts >= min_count).sum(axis=0) >= smallest_group
    elif policy == PREFILTER_ALL_REPLICATES_ONE_CONDITION:
        keep = counts.columns.to_series(index=counts.columns).map(lambda _: False)
        for condition in group_sizes.index:
            members = conditions.index[conditions.eq(condition)]
            keep |= counts.loc[members].ge(min_count).all(axis=0)
    else:
        raise ValidationError(
            f"Unknown DA prefilter policy {policy!r}; expected one of {', '.join(PREFILTER_POLICIES)}."
        )
    return keep, smallest_group


def _matrix_to_analysis_counts(raw_matrix, sample_ids: list[str], sample_to_analysis: dict[str, str], pd):
    """Transpose a matrix and sum declared technical replicates into model rows."""

    absent = [sample_id for sample_id in sample_ids if sample_id not in raw_matrix.columns]
    if absent:
        raise ValidationError(
            f"Count matrix is missing target samples: {', '.join(absent)}."
        )
    per_library = raw_matrix.set_index("unit_id")[sample_ids].T
    per_library.index.name = "sample_id"
    per_library["analysis_id"] = [sample_to_analysis[sample_id] for sample_id in per_library.index]
    return per_library.groupby("analysis_id", sort=True).sum(numeric_only=True)


def read_library_depths(path: Path, expected_sample_ids: Iterable[str]) -> dict[str, int]:
    """Read audited retained-library depths for fixed external normalization."""

    if not path.is_file():
        raise ValidationError(f"Library summary does not exist: {path}")
    expected = set(expected_sample_ids)
    depths: dict[str, int] = {}
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        depth_column = (
            "normalization_depth"
            if "normalization_depth" in fields
            else "input_fastq_reads"
        )
        required = {"sample_id", depth_column}
        missing = sorted(required - fields)
        if missing:
            raise ValidationError(
                f"{path}: library summary is missing columns: {', '.join(missing)}."
            )
        for row_number, row in enumerate(reader, start=2):
            sample_id = (row.get("sample_id") or "").strip()
            if not sample_id:
                raise ValidationError(f"{path}:{row_number}: blank sample_id.")
            if sample_id in depths:
                raise ValidationError(f"{path}:{row_number}: duplicate sample_id {sample_id!r}.")
            try:
                depth = int((row.get(depth_column) or "").strip())
            except ValueError as exc:
                raise ValidationError(
                    f"{path}:{row_number}: {depth_column} must be a positive integer."
                ) from exc
            if depth <= 0:
                raise ValidationError(
                    f"{path}:{row_number}: {depth_column} must be a positive integer."
                )
            depths[sample_id] = depth
    absent = sorted(expected - set(depths))
    unexpected = sorted(set(depths) - expected)
    if absent or unexpected:
        details = []
        if absent:
            details.append(f"missing target samples: {', '.join(absent)}")
        if unexpected:
            details.append(f"unexpected samples: {', '.join(unexpected)}")
        raise ValidationError(f"{path}: library-depth sample mismatch ({'; '.join(details)}).")
    return depths


def analysis_depths_and_size_factors(
    *,
    sample_to_analysis: dict[str, str],
    library_depths: dict[str, int],
) -> tuple[dict[str, int], dict[str, float]]:
    """Combine technical-library depths and scale fixed factors to geometric mean one."""

    analysis_depths: dict[str, int] = defaultdict(int)
    for sample_id, analysis_id in sample_to_analysis.items():
        try:
            analysis_depths[analysis_id] += library_depths[sample_id]
        except KeyError as exc:
            raise ValidationError(f"No audited library depth is available for {sample_id!r}.") from exc
    geometric_mean = math.exp(
        sum(math.log(depth) for depth in analysis_depths.values()) / len(analysis_depths)
    )
    size_factors = {
        analysis_id: depth / geometric_mean
        for analysis_id, depth in analysis_depths.items()
    }
    return dict(analysis_depths), size_factors


def _run_deseq2_with_fixed_size_factors(dds, size_factors, np) -> str:
    """Run the pinned PyDESeq2 fit and return the dispersion trend actually used."""

    ordered = np.asarray([size_factors[str(index)] for index in dds.obs_names], dtype=float)
    dds.obs["size_factors"] = ordered
    normed_counts = np.asarray(dds.X, dtype=float) / ordered[:, None]
    dds.layers["normed_counts"] = normed_counts
    dds.var["_normed_means"] = normed_counts.mean(axis=0)
    dds.fit_genewise_dispersions()
    dds.fit_dispersion_trend()
    dds.fit_dispersion_prior()
    dds.fit_MAP_dispersions()
    dds.fit_LFC()
    dds.calculate_cooks()
    if dds.refit_cooks:
        dds.refit()
    dds.cooks_outlier()
    return str(dds.uns.get("disp_function_type", "unknown"))


def _write_skip(outdir: Path, reason: str, metadata: dict[str, object]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "SKIPPED.txt").write_text(reason + "\n", encoding="utf-8")
    with (outdir / "model_metadata.json").open("wt", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _run_one_matrix(
    path: Path,
    *,
    samples: tuple[Sample, ...],
    outdir: Path,
    min_count: int,
    alpha: float,
    threads: int,
    library_depths: dict[str, int],
    library_summary_path: Path,
    contrasts: tuple[tuple[str, str], ...],
    effective_analysis_depths: dict[str, float] | None = None,
    normalization_method: str | None = None,
    normalization_source: Path | None = None,
    prefilter_policy: str = PREFILTER_ANY_SAMPLES,
) -> None:
    np, pd, DeseqDataSet, DeseqStats = _optional_dependencies()
    raw = pd.read_csv(path, sep="\t", dtype={"unit_id": str})
    if "unit_id" not in raw.columns:
        raise ValidationError(f"{path}: no unit_id column.")
    sample_ids = [sample.sample_id for sample in samples]
    sample_to_analysis, metadata = _combine_technical_replicates(samples, pd)
    counts = _matrix_to_analysis_counts(raw, sample_ids, sample_to_analysis, pd)
    metadata = metadata.loc[counts.index]
    analysis_depths, size_factors = analysis_depths_and_size_factors(
        sample_to_analysis=sample_to_analysis,
        library_depths=library_depths,
    )
    analysis_depths = {analysis_id: analysis_depths[analysis_id] for analysis_id in counts.index}
    size_factors = {analysis_id: size_factors[analysis_id] for analysis_id in counts.index}
    if effective_analysis_depths is not None:
        if not normalization_method or normalization_source is None or not normalization_source.is_file():
            raise ValidationError("External normalization requires a method and an existing provenance file.")
        size_factors = external_size_factors(effective_analysis_depths, analysis_depths)
    elif normalization_method is not None or normalization_source is not None:
        raise ValidationError("External normalization metadata requires effective analysis depths.")
    conditions = metadata["condition"].astype(str)
    condition_sizes = conditions.value_counts().sort_index()
    model_info: dict[str, object] = {
        "input_count_matrix": str(path),
        "library_summary": str(library_summary_path),
        "size_factor_method": normalization_method or "fixed_total_retained_library_abundance",
        "dispersion_fit_type_requested": "parametric",
        "normalization_scope": "one shared sample-specific factor across all count projections",
        "reporting_normalization": f"CP{REPORTING_DEPTH // 1_000_000}M",
        "analysis_library_depths": analysis_depths,
        "size_factors_geometric_mean_one": size_factors,
        "fdr_scope": "within projection and contrast",
        "alpha": alpha,
        "min_count": min_count,
        "prefilter_policy": prefilter_policy,
        "condition_sizes": {key: int(value) for key, value in condition_sizes.items()},
        "n_technical_libraries_combined": int(len(samples) - len(counts)),
        "declared_contrasts": [
            {"numerator": numerator, "denominator": denominator}
            for numerator, denominator in contrasts
        ],
    }
    if effective_analysis_depths is not None:
        model_info.update(
            effective_analysis_depths=effective_analysis_depths,
            normalization_source=str(normalization_source),
            reporting_normalization="CP30M uses original retained depth; effective_CP30M uses external effective depth",
        )
    if len(condition_sizes) < 2:
        _write_skip(outdir, "Only one condition is present; no pairwise DA contrast is possible.", model_info)
        return
    if int(condition_sizes.min()) < 2:
        _write_skip(
            outdir,
            "At least one condition has fewer than two independent biological replicates; DA is withheld.",
            model_info,
        )
        return
    keep, smallest_group = _filter_features(counts, conditions, min_count, prefilter_policy)
    model_info["smallest_group_size"] = smallest_group
    model_info["n_features_before_prefilter"] = int(counts.shape[1])
    model_info["n_features_after_prefilter"] = int(keep.sum())
    if prefilter_policy == PREFILTER_ALL_REPLICATES_ONE_CONDITION:
        rejected_reason = f"not >= {min_count} raw counts in every replicate of any one condition"
    else:
        rejected_reason = f"fewer than {smallest_group} samples with >= {min_count} raw counts"
    feature_status = pd.DataFrame(
        {
            "unit_id": counts.columns,
            "prefilter_keep": keep.to_numpy(),
            "prefilter_reason": [
                "kept" if keep.loc[unit_id] else rejected_reason
                for unit_id in counts.columns
            ],
        }
    )
    outdir.mkdir(parents=True, exist_ok=True)
    feature_status.to_csv(outdir / "feature_filtering.tsv", sep="\t", index=False)
    if int(keep.sum()) < MIN_FEATURES_FOR_DA:
        _write_skip(
            outdir,
            f"Only {int(keep.sum())} features pass prefilter; need at least {MIN_FEATURES_FOR_DA} for a stable dispersion fit.",
            model_info,
        )
        return
    filtered_counts = counts.loc[:, keep]
    design, design_warning = _choose_design(metadata, np, pd)
    model_info["design"] = design
    if design_warning:
        model_info["design_warning"] = design_warning
    dds = DeseqDataSet(
        counts=filtered_counts.astype(int),
        metadata=metadata,
        design=design,
        refit_cooks=True,
        n_cpus=threads,
        quiet=True,
    )
    model_info["dispersion_fit_type_used"] = _run_deseq2_with_fixed_size_factors(
        dds, size_factors, np
    )
    model_info["dispersion_fallback_used"] = (
        model_info["dispersion_fit_type_used"] != model_info["dispersion_fit_type_requested"]
    )
    normalized = None
    if "normed_counts" in dds.layers:
        normalized = pd.DataFrame(
            dds.layers["normed_counts"], index=filtered_counts.index, columns=filtered_counts.columns
        )
    elif "size_factors" in dds.obsm:
        normalized = filtered_counts.div(dds.obsm["size_factors"], axis=0)
    if normalized is not None:
        normalized.to_csv(outdir / "normalized_counts.tsv", sep="\t", index_label="analysis_id")
    depth_series = pd.Series(analysis_depths, dtype=float).loc[filtered_counts.index]
    cp30m = filtered_counts.div(depth_series, axis=0) * REPORTING_DEPTH
    cp30m.to_csv(outdir / "cp30m_counts.tsv", sep="\t", index_label="analysis_id")
    if effective_analysis_depths is not None:
        effective_depth_series = pd.Series(effective_analysis_depths).loc[filtered_counts.index]
        effective_cp30m = filtered_counts.div(effective_depth_series, axis=0) * REPORTING_DEPTH
        effective_cp30m.to_csv(outdir / "effective_cp30m_counts.tsv", sep="\t", index_label="analysis_id")
    pd.DataFrame(
        {
            "analysis_id": filtered_counts.index,
            "normalization_depth": [analysis_depths[index] for index in filtered_counts.index],
            "deseq2_size_factor": [size_factors[index] for index in filtered_counts.index],
            "reporting_scale": [f"CP{REPORTING_DEPTH // 1_000_000}M"] * len(filtered_counts.index),
        }
    ).to_csv(outdir / "normalization_factors.tsv", sep="\t", index=False)
    metadata.to_csv(outdir / "analysis_sample_metadata.tsv", sep="\t", index_label="analysis_id")
    for numerator, denominator in contrasts:
        stats = DeseqStats(
            dds,
            contrast=["condition", numerator, denominator],
            alpha=alpha,
            independent_filter=True,
            quiet=True,
            n_cpus=threads,
        )
        stats.summary()
        result = stats.results_df.copy()
        result.index.name = "unit_id"
        result = result.reset_index()
        feature_metadata = raw.drop(columns=sample_ids).drop_duplicates("unit_id")
        result = result.merge(feature_metadata, on="unit_id", how="left", validate="one_to_one")
        result.to_csv(
            outdir / f"{numerator}_vs_{denominator}.tsv", sep="\t", index=False
        )
    with (outdir / "model_metadata.json").open("wt", encoding="utf-8") as handle:
        json.dump(model_info, handle, indent=2, sort_keys=True)
        handle.write("\n")


def external_size_factors(effective_depths: dict[str, float], expected_analysis_ids) -> dict[str, float]:
    """Validate externally estimated effective depths; never replace integer raw counts."""
    expected = set(expected_analysis_ids)
    if not expected or set(effective_depths) != expected:
        raise ValidationError("External normalization must contain exactly the analysis sample IDs.")
    if any(not math.isfinite(float(value)) or float(value) <= 0 for value in effective_depths.values()):
        raise ValidationError("External effective depths must be positive and finite.")
    center = math.exp(sum(math.log(value) for value in effective_depths.values()) / len(expected))
    return {key: float(effective_depths[key]) / center for key in expected_analysis_ids}


def run_differential_abundance(
    *,
    count_directory: Path,
    samples: tuple[Sample, ...],
    outdir: Path,
    min_count: int,
    alpha: float,
    threads: int,
    library_summary_path: Path,
    contrasts: tuple[tuple[str, str], ...],
    prefilter_policy: str = PREFILTER_ANY_SAMPLES,
) -> None:
    """Run a separately labelled PyDESeq2 model for every raw count projection.

    Matrices remain integer-valued because ambiguous reads are represented by
    named equivalence groups rather than fractional allocation. Each projection
    is an alternate analysis resolution and therefore receives its own model;
    this function never treats locus and product results as independent
    confirmation of each other.
    """

    matrix_paths = sorted(count_directory.glob("*.raw.tsv"))
    if not matrix_paths:
        raise ValidationError(f"No raw count matrices found under {count_directory}.")
    sample_ids = [sample.sample_id for sample in samples]
    library_depths = read_library_depths(library_summary_path, sample_ids)
    for matrix_path in matrix_paths:
        projection = matrix_path.name[: -len(".raw.tsv")]
        print(f"[{projection}] DA: {len(contrasts)} declared contrasts", file=sys.stderr, flush=True)
        _run_one_matrix(
            matrix_path,
            samples=samples,
            outdir=outdir / projection,
            min_count=min_count,
            alpha=alpha,
            threads=threads,
            library_depths=library_depths,
            library_summary_path=library_summary_path,
            contrasts=contrasts,
            prefilter_policy=prefilter_policy,
        )
        print(f"[{projection}] DA finished (see model metadata for fit or skip status)", file=sys.stderr, flush=True)
