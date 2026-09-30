from __future__ import annotations

import gzip
import json
from argparse import Namespace

import pytest

from phasis_srna_da.cli import (
    _validate_declared_contrasts,
    _validate_run_policy,
    build_parser,
    main,
)
from phasis_srna_da.errors import ValidationError


def _shared_arguments() -> list[str]:
    return [
        "--processed-libraries",
        "libraries",
        "--library-format",
        "fastq",
        "--min-length",
        "18",
        "--max-length",
        "30",
        "--targets",
        "targets.tsv",
        "--assembly-id",
        "Zm-B73-REFERENCE-NAM-5.0",
        "--annotation-id",
        "frozen-sRNA-test-v1",
        "--reference",
        "miRNA=mature_mirna.fa",
        "--phasiRNAs",
        "21_phasiRNAs.tsv",
        "--reference-orientation",
        "molecule",
    ]


def test_cli_reports_development_version(capsys):
    with pytest.raises(SystemExit, match="0"):
        build_parser().parse_args(["--version"])

    assert capsys.readouterr().out.strip() == "phasis-srna-da 0.2.0.dev0"


def test_cli_requires_phasis_products_and_explicit_reference_orientation():
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "validate",
                "--processed-libraries",
                "libraries",
                "--library-format",
                "fastq",
                "--min-length",
                "18",
                "--max-length",
                "30",
                "--targets",
                "targets.tsv",
                "--assembly-id",
                "Zm-B73-REFERENCE-NAM-5.0",
                "--annotation-id",
                "frozen-sRNA-test-v1",
                "--reference",
                "miRNA=mature_mirna.fa",
                "--reference-orientation",
                "molecule",
            ]
        )
    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "validate",
                "--processed-libraries",
                "libraries",
                "--library-format",
                "fastq",
                "--min-length",
                "18",
                "--max-length",
                "30",
                "--targets",
                "targets.tsv",
                "--assembly-id",
                "Zm-B73-REFERENCE-NAM-5.0",
                "--annotation-id",
                "frozen-sRNA-test-v1",
                "--reference",
                "miRNA=mature_mirna.fa",
                "--phasiRNAs",
                "21_phasiRNAs.tsv",
            ]
        )

    args = parser.parse_args(["validate", *_shared_arguments()])
    assert [path.name for path in args.phasi_catalog] == ["21_phasiRNAs.tsv"]
    assert args.reference_orientation == "molecule"


def test_legacy_phasi_catalog_spelling_remains_an_alias():
    arguments = _shared_arguments()
    option_index = arguments.index("--phasiRNAs")
    arguments[option_index] = "--phasi-catalog"

    args = build_parser().parse_args(["validate", *arguments])

    assert [path.name for path in args.phasi_catalog] == ["21_phasiRNAs.tsv"]


def test_da_requires_an_explicit_supported_normalization_policy():
    with pytest.raises(ValidationError, match="--normalization is required"):
        _validate_run_policy(Namespace(run_da=True, normalization=None, contrast=[]))

    with pytest.raises(ValidationError, match="only valid with --run-da"):
        _validate_run_policy(
            Namespace(run_da=False, normalization="total-qc-depth", contrast=[])
        )

    with pytest.raises(ValidationError, match="--contrast is required"):
        _validate_run_policy(
            Namespace(run_da=True, normalization="total-qc-depth", contrast=[])
        )

    _validate_run_policy(
        Namespace(
            run_da=True,
            normalization="total-qc-depth",
            contrast=[("ms28_2", "WT_2")],
        )
    )


def test_cli_parses_repeated_direct_contrasts_and_validates_conditions():
    args = build_parser().parse_args(
        [
            "run",
            *_shared_arguments(),
            "--outdir",
            "run",
            "--run-da",
            "--normalization",
            "total-qc-depth",
            "--contrast",
            "ms28_0p4,WT_0p4",
            "--contrast",
            "ms28_2,WT_2",
        ]
    )

    assert args.contrast == [("ms28_0p4", "WT_0p4"), ("ms28_2", "WT_2")]
    _validate_declared_contrasts(
        args.contrast, {"ms28_0p4", "WT_0p4", "ms28_2", "WT_2"}
    )
    with pytest.raises(ValidationError, match="not present in --targets"):
        _validate_declared_contrasts(args.contrast, {"ms28_0p4", "WT_0p4"})


def test_validate_reports_required_phasis_and_orientation_contract(tmp_path, capsys):
    libraries = tmp_path / "libraries"
    libraries.mkdir()
    with gzip.open(libraries / "sample_1.fastq.gz", "wt", encoding="utf-8") as handle:
        handle.write("@read_1\nAAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIIII\n")
    targets = tmp_path / "targets.tsv"
    targets.write_text(
        "sample_id\tcondition\tbiological_replicate\n"
        "sample_1\tcontrol\t1\n",
        encoding="utf-8",
    )
    reference = tmp_path / "mature_mirna.fa"
    reference.write_text(">miR-test\nCCCCCCCCCCCCCCCCCCCCC\n", encoding="utf-8")
    phasi_products = tmp_path / "21_phasiRNAs.tsv"
    phasi_products.write_text(
        "identifier\tphase\twindow_unit_id\tstrand\texpected_register_pos\t"
        "register_class\ttag_seq\n"
        "1:1..100\t21\tmain\tw\t21\tcore_exact\tAAAAAAAAAAAAAAAAAAAAA\n",
        encoding="utf-8",
    )

    exit_code = main(
        [
            "validate",
            "--processed-libraries",
            str(libraries),
            "--library-format",
            "fastq",
            "--min-length",
            "18",
            "--max-length",
            "30",
            "--targets",
            str(targets),
            "--assembly-id",
            "Zm-B73-REFERENCE-NAM-5.0",
            "--annotation-id",
            "frozen-sRNA-test-v1",
            "--reference",
            f"miRNA={reference}",
            "--phasiRNAs",
            str(phasi_products),
            "--reference-orientation",
            "molecule",
            "--json",
        ]
    )

    summary = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert summary["phasi_product_files"] == 1
    assert summary["phasi_catalog_rows_selected"] == 1
    assert summary["assembly_id"] == "Zm-B73-REFERENCE-NAM-5.0"
    assert summary["annotation_id"] == "frozen-sRNA-test-v1"
    assert summary["reference_orientation"] == "molecule"
    assert summary["mapping_orientation"] == "direct_only_no_reverse_complement"


def test_tiny_tag_run_drops_n_rows_preserves_phasis_outputs_and_records_contract(tmp_path):
    libraries = tmp_path / "libraries"
    libraries.mkdir()
    (libraries / "sample_1.tag").write_text(
        "AAAAAAAAAAAAAAAAAAAAA\t5\n"
        "AAAAANAAAAAAAAAAAAAAA\t7\n",
        encoding="utf-8",
    )
    targets = tmp_path / "targets.tsv"
    targets.write_text(
        "sample_id\tcondition\tbiological_replicate\n"
        "sample_1\tcontrol\t1\n",
        encoding="utf-8",
    )
    reference = tmp_path / "mature_mirna.fa"
    reference.write_text(">miR-test\nCCCCCCCCCCCCCCCCCCCCC\n", encoding="utf-8")
    phasi_products = tmp_path / "21_phasiRNAs.tsv"
    phasi_products.write_text(
        "identifier\tphase\twindow_unit_id\tstrand\texpected_register_pos\t"
        "register_class\ttag_seq\n"
        "1:1..100\t21\tmain\tw\t21\tcore_exact\tAAAAAAAAAAAAAAAAAAAAA\n",
        encoding="utf-8",
    )
    fake_bowtie_build = tmp_path / "bowtie-build"
    fake_bowtie_build.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "if '--version' in sys.argv:\n"
        "    print('bowtie-build version 1-test')\n",
        encoding="utf-8",
    )
    fake_bowtie_build.chmod(0o755)
    fake_bowtie = tmp_path / "bowtie"
    fake_bowtie.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "if '--version' in sys.argv:\n"
        "    print('bowtie version 1-test')\n"
        "else:\n"
        "    print('@HD\\tVN:1.0')\n"
        "    print('read_1\\t0\\tref_000000002\\t1\\t255\\t21M\\t*\\t0\\t0\\t'"
        " + 'AAAAAAAAAAAAAAAAAAAAA' + '\\t*')\n",
        encoding="utf-8",
    )
    fake_bowtie.chmod(0o755)
    outdir = tmp_path / "run"

    exit_code = main(
        [
            "run",
            "--processed-libraries",
            str(libraries),
            "--library-format",
            "tag",
            "--min-length",
            "18",
            "--max-length",
            "30",
            "--targets",
            str(targets),
            "--assembly-id",
            "Zm-B73-REFERENCE-NAM-5.0",
            "--annotation-id",
            "frozen-sRNA-test-v1",
            "--reference",
            f"miRNA={reference}",
            "--phasiRNAs",
            str(phasi_products),
            "--reference-orientation",
            "molecule",
            "--bowtie",
            str(fake_bowtie),
            "--bowtie-build",
            str(fake_bowtie_build),
            "--outdir",
            str(outdir),
        ]
    )

    metadata = json.loads(
        (outdir / "provenance" / "run_metadata.json").read_text(encoding="utf-8")
    )
    assert exit_code == 0
    assert (outdir / "counts" / "phasi_product.raw.tsv").is_file()
    assert (outdir / "counts" / "phasi_locus.raw.tsv").is_file()
    assert "ref_000000002" in (
        outdir / "counts" / "phasi_product.raw.tsv"
    ).read_text(encoding="utf-8")
    assert metadata["status"] == "complete"
    assert metadata["assembly_id"] == "Zm-B73-REFERENCE-NAM-5.0"
    assert metadata["annotation_id"] == "frozen-sRNA-test-v1"
    assert metadata["reference_orientation"] == "molecule"
    assert metadata["mapping_orientation"] == "direct_only_no_reverse_complement"
    assert metadata["library_format"] == "tag"
    assert metadata["ambiguous_N_policy"] == "drop_complete_library_record"
    assert metadata["da_normalization_policy"] is None
    assert metadata["command_line"][:2] == ["phasis-srna-da", "run"]
    assert str(phasi_products) in metadata["command_line"]
    assert f"miRNA={reference}" in metadata["command_line"]
    summary = (outdir / "qc" / "library_summary.tsv").read_text(encoding="utf-8")
    assert "excluded_N_records\texcluded_N_abundance" in summary
    assert "\t2\t12\t1\t7\t0\t0\t1\t1\t5\t5\t" in summary
