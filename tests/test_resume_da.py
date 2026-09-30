from __future__ import annotations

import json

import pytest

from phasis_srna_da.cli import main
from phasis_srna_da.errors import ValidationError
from phasis_srna_da.io import sha256_file, write_tsv
from phasis_srna_da.provenance import write_json
from phasis_srna_da.resume_da import REQUIRED_PROJECTIONS, check_count_run


@pytest.fixture
def count_run(tmp_path):
    source = tmp_path / "counts_run"
    write_json(source / "provenance/run_metadata.json", {
        "tool": "phasis-srna-da", "status": "complete", "assembly_id": "test-assembly",
        "annotation_id": "frozen-test", "library_format": "tag", "library_min_length": 31,
        "library_max_length": 50, "ambiguous_N_policy": "drop_complete_library_record",
        "reference_orientation": "molecule", "mapping_orientation": "direct_only_no_reverse_complement",
    })
    write_tsv(source / "provenance/input_manifest.tsv", ["role", "path", "sha256"], [
        {"role": "phasis_product_catalog", "path": "/original/21_phasiRNAs.tsv", "sha256": "a" * 64}
    ])
    samples = [
        {"sample_id": f"{condition}_{rep}", "condition": condition, "biological_replicate": str(rep)}
        for condition in ("WT", "ms28") for rep in (1, 2, 3)
    ]
    write_tsv(source / "provenance/resolved_targets.tsv", list(samples[0]), samples)
    qc = [{**row, "retained_abundance": 100, "normalization_depth": 100,
           "primary_unique_assigned_abundance": 20, "primary_equivalence_assigned_abundance": 0,
           "primary_unassigned_abundance": 80} for row in samples]
    write_tsv(source / "qc/library_summary.tsv", list(qc[0]), qc)
    conservation = []
    for projection in sorted(REQUIRED_PROJECTIONS):
        assigned = 20 if projection == "global_atomic" else 0
        rows = [{"unit_id": "tRNA_unit", "assignment_type": "unique", **{s["sample_id"]: assigned for s in samples}}] if assigned else []
        write_tsv(source / f"counts/{projection}.raw.tsv", ["unit_id", "assignment_type", *[s["sample_id"] for s in samples]], rows)
        conservation.extend({
            "projection": projection, "sample_id": s["sample_id"], "input_abundance": 100,
            "unique_assigned_abundance": assigned, "equivalence_assigned_abundance": 0,
            "assigned_abundance": assigned, "unassigned_abundance": 100 - assigned, "conserved": "true",
        } for s in samples)
    write_tsv(source / "mapping/conservation_by_sample.tsv", list(conservation[0]), conservation)
    return source


def test_completed_count_run_accepts_empty_phasi_matrices_for_long_read_stratum(count_run):
    metadata, samples, paths = check_count_run(count_run)
    assert metadata["library_min_length"] == 31
    assert len(samples) == 6
    assert len(paths) == 10


@pytest.mark.parametrize("change,match", [
    ("running", "completed"), ("missing_phasis", "mandatory Phasis"),
    ("matrix_count", "conservation"), ("negative_count", "nonnegative integers"),
    ("fractional_count", "nonnegative integers"), ("missing_matrix", "mandatory global/Phasis"),
    ("depth", "normalization depth"), ("duplicate_unit", "duplicate count unit"),
    ("condition", "design disagree"),
])
def test_bad_count_evidence_is_rejected(count_run, change, match):
    matrix = count_run / "counts/global_atomic.raw.tsv"
    if change == "running":
        path = count_run / "provenance/run_metadata.json"
        meta = json.loads(path.read_text())
        meta["status"] = "running"
        write_json(path, meta)
    elif change == "missing_phasis":
        (count_run / "provenance/input_manifest.tsv").write_text("role\tpath\tsha256\n")
    elif change in {"matrix_count", "negative_count", "fractional_count"}:
        replacement = {"matrix_count": "21", "negative_count": "-1", "fractional_count": "20.5"}[change]
        matrix.write_text(matrix.read_text().replace("\t20", f"\t{replacement}", 1))
    elif change == "missing_matrix":
        (count_run / "counts/phasi_product.raw.tsv").unlink()
    elif change == "depth":
        path = count_run / "qc/library_summary.tsv"
        path.write_text(path.read_text().replace("\t100\t100", "\t100\t101", 1))
    elif change == "condition":
        path = count_run / "qc/library_summary.tsv"
        path.write_text(path.read_text().replace("WT_1\tWT", "WT_1\tms28", 1))
    elif change == "duplicate_unit":
        text = matrix.read_text()
        matrix.write_text(text + text.splitlines()[1] + "\n")
    with pytest.raises(ValidationError, match=match):
        check_count_run(count_run)


def test_da_reuses_counts_preserves_source_and_passes_explicit_contrasts(count_run, tmp_path, monkeypatch):
    called = []
    monkeypatch.setattr("phasis_srna_da.da.run_differential_abundance", lambda **kwargs: called.append(kwargs))
    monkeypatch.setattr("phasis_srna_da.da._optional_dependencies", lambda: None)
    monkeypatch.setattr("phasis_srna_da.resume_da.version", lambda name: "test-version")
    before = {p: sha256_file(p) for p in count_run.rglob("*") if p.is_file()}
    output = tmp_path / "da_run"
    assert main(["da", "--count-run", str(count_run), "--outdir", str(output),
                 "--normalization", "total-qc-depth", "--contrast", "ms28,WT"]) == 0
    assert len(called) == 1
    assert called[0]["contrasts"] == (("ms28", "WT"),)
    assert called[0]["library_summary_path"] == count_run / "qc/library_summary.tsv"
    assert called[0]["count_directory"] == count_run / "counts"
    assert called[0]["prefilter_policy"] == "all-replicates-in-one-condition"
    assert before == {p: sha256_file(p) for p in count_run.rglob("*") if p.is_file()}
    metadata = json.loads((output / "provenance/run_metadata.json").read_text())
    assert metadata["status"] == "complete"
    assert metadata["library_min_length"] == 31
    assert metadata["annotation_id"] == "frozen-test"
    assert metadata["da_contrasts"] == [{"numerator": "ms28", "denominator": "WT"}]
    assert metadata["prefilter_policy"] == "all-replicates-in-one-condition"
    assert metadata["command_line"][:2] == ["phasis-srna-da", "da"]
    assert "source_count_evidence" in (output / "provenance/input_manifest.tsv").read_text()


@pytest.mark.parametrize("bad_args", [
    ["--contrast", "missing,WT"], ["--contrast", "ms28,WT", "--contrast", "ms28,WT"],
    ["--contrast", "ms28,WT", "--contrast", "WT,ms28"], ["--contrast", "ms28,WT", "--fdr", "2"],
])
def test_da_rejects_bad_contrasts_before_creating_outputs(count_run, tmp_path, bad_args):
    output = tmp_path / "bad_da"
    with pytest.raises(SystemExit, match="2"):
        main(["da", "--count-run", str(count_run), "--outdir", str(output),
              "--normalization", "total-qc-depth", *bad_args])
    assert not output.exists()


def test_da_rejects_output_inside_count_run(count_run):
    output = count_run / "new_da"
    with pytest.raises(SystemExit, match="2"):
        main(["da", "--count-run", str(count_run), "--outdir", str(output),
              "--normalization", "total-qc-depth", "--contrast", "ms28,WT"])
    assert not output.exists()


def test_da_failure_is_recorded_and_never_overwritten(count_run, tmp_path, monkeypatch):
    def fail(**kwargs):
        raise ValidationError("test fit failure")
    monkeypatch.setattr("phasis_srna_da.da.run_differential_abundance", fail)
    monkeypatch.setattr("phasis_srna_da.da._optional_dependencies", lambda: None)
    monkeypatch.setattr("phasis_srna_da.resume_da.version", lambda name: "test-version")
    output = tmp_path / "failed_da"
    args = ["da", "--count-run", str(count_run), "--outdir", str(output),
            "--normalization", "total-qc-depth", "--contrast", "ms28,WT"]
    with pytest.raises(SystemExit, match="2"):
        main(args)
    metadata_path = output / "provenance/run_metadata.json"
    assert json.loads(metadata_path.read_text())["status"] == "failed"
    before = sha256_file(metadata_path)
    with pytest.raises(SystemExit, match="2"):
        main(args)
    assert sha256_file(metadata_path) == before
