from __future__ import annotations

import math
import json

import pytest

from phasis_srna_da.da import (
    PREFILTER_ALL_REPLICATES_ONE_CONDITION,
    PREFILTER_ANY_SAMPLES,
    _filter_features,
    analysis_depths_and_size_factors,
    read_library_depths,
)
from phasis_srna_da.errors import ValidationError
from phasis_srna_da.io import write_tsv
from phasis_srna_da.models import Sample


def test_reads_audited_retained_library_depths_and_rejects_sample_mismatch(tmp_path):
    summary = tmp_path / "library_summary.tsv"
    summary.write_text(
        "sample_id\tnormalization_depth\nA\t100\nB\t400\n",
        encoding="utf-8",
    )

    assert read_library_depths(summary, ["A", "B"]) == {"A": 100, "B": 400}
    with pytest.raises(ValidationError, match="sample mismatch"):
        read_library_depths(summary, ["A", "C"])


def test_reads_legacy_fastq_depth_column_for_compatible_old_results(tmp_path):
    summary = tmp_path / "library_summary.tsv"
    summary.write_text(
        "sample_id\tinput_fastq_reads\nA\t100\nB\t400\n",
        encoding="utf-8",
    )

    assert read_library_depths(summary, ["A", "B"]) == {"A": 100, "B": 400}


def test_combines_technical_depths_and_scales_size_factors_to_geometric_mean_one():
    depths, factors = analysis_depths_and_size_factors(
        sample_to_analysis={"A_run1": "A", "A_run2": "A", "B": "B"},
        library_depths={"A_run1": 40, "A_run2": 60, "B": 400},
    )

    assert depths == {"A": 100, "B": 400}
    assert factors == pytest.approx({"A": 0.5, "B": 2.0})
    assert math.prod(factors.values()) ** (1 / len(factors)) == pytest.approx(1.0)


def test_condition_replicate_prefilter_requires_one_complete_biological_group():
    pd = pytest.importorskip("pandas")
    index = [f"{condition}_{rep}" for condition in ("A", "B") for rep in (1, 2, 3)]
    conditions = pd.Series(["A"] * 3 + ["B"] * 3, index=index)
    counts = pd.DataFrame(
        {
            "scattered": [10, 0, 0, 10, 10, 0],
            "complete_A": [10, 11, 12, 0, 0, 0],
            "complete_B": [0, 0, 0, 10, 10, 10],
            "below": [9, 9, 9, 10, 10, 9],
        },
        index=index,
    )
    any_samples, _ = _filter_features(counts, conditions, 10, PREFILTER_ANY_SAMPLES)
    one_condition, _ = _filter_features(
        counts, conditions, 10, PREFILTER_ALL_REPLICATES_ONE_CONDITION
    )
    assert any_samples.to_dict() == {
        "scattered": True, "complete_A": True, "complete_B": True, "below": False
    }
    assert one_condition.to_dict() == {
        "scattered": False, "complete_A": True, "complete_B": True, "below": False
    }


def test_real_da_fits_only_nine_declared_ms28_contrasts(tmp_path):
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")
    pytest.importorskip("pydeseq2")
    from phasis_srna_da.da import run_differential_abundance

    conditions = [f"{genotype}_{stage}" for genotype in ("WT", "ms28") for stage in ("0p4", "2", "5")]
    samples = tuple(Sample(sample_id=f"{c}_{r}", condition=c, biological_replicate=str(r),
                           library_id=f"{c}_{r}") for c in conditions for r in (1, 2, 3))
    contrasts = tuple([(f"ms28_{s}", f"WT_{s}") for s in ("0p4", "2", "5")] +
                      [(f"{g}_{n}", f"{g}_{d}") for g in ("WT", "ms28")
                       for n, d in (("2", "0p4"), ("5", "0p4"), ("5", "2"))])
    rng = np.random.default_rng(741)
    counts = rng.negative_binomial(20, 0.08, size=(40, len(samples)))
    rows = [{"unit_id": f"unit_{i}", **{s.sample_id: int(counts[i, j]) for j, s in enumerate(samples)}}
            for i in range(len(counts))]
    write_tsv(tmp_path / "counts/global_atomic.raw.tsv", list(rows[0]), rows)
    depths = [{"sample_id": s.sample_id, "normalization_depth": (i + 1) * 100000}
              for i, s in enumerate(samples)]
    summary = tmp_path / "library_summary.tsv"
    write_tsv(summary, list(depths[0]), depths)
    run_differential_abundance(count_directory=tmp_path / "counts", samples=samples,
                              outdir=tmp_path / "da", min_count=10, alpha=0.05,
                              threads=1, library_summary_path=summary, contrasts=contrasts)
    output = tmp_path / "da/global_atomic"
    assert {p.name for p in output.glob("*_vs_*.tsv")} == {f"{n}_vs_{d}.tsv" for n, d in contrasts}
    assert len(contrasts) == 9
    metadata = json.loads((output / "model_metadata.json").read_text())
    assert metadata["size_factor_method"] == "fixed_total_retained_library_abundance"
    assert len(metadata["declared_contrasts"]) == 9
    assert metadata["condition_sizes"] == dict.fromkeys(conditions, 3)
    factors = pd.read_csv(output / "normalization_factors.tsv", sep="\t")
    assert set(factors["normalization_depth"]) == {r["normalization_depth"] for r in depths}
