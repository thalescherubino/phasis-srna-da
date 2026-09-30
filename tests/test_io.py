from __future__ import annotations

import gzip

import pytest

from phasis_srna_da.errors import ValidationError
from phasis_srna_da.io import (
    find_processed_libraries,
    read_collapsed_tag,
    read_collapsed_fastq,
    read_fastq,
    read_fastq_collapsed,
    read_targets,
    resolve_target_libraries,
    write_mapping_queries,
)


VALID_LIBRARY = (
    "@read_one\nACGTACGTACGTACGTACGTA\n+\nIIIIIIIIIIIIIIIIIIIII\n"
    "@read_two\nTGCATGCATGCATGCATGCAT\n+\nIIIIIIIIIIIIIIIIIIIII\n"
)


def test_reads_quality_controlled_fastq_and_assigns_one_count_per_read(tmp_path):
    path = tmp_path / "sample.fastq.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(VALID_LIBRARY)

    records = read_fastq(path)

    assert [(record.query_id, record.sequence) for record in records] == [
        ("read_1", "ACGTACGTACGTACGTACGTA"),
        ("read_2", "TGCATGCATGCATGCATGCAT"),
    ]
    queries = tmp_path / "queries.fasta"
    write_mapping_queries(records, queries)
    assert queries.read_text(encoding="utf-8").startswith(">read_1\nACGTACGTACGTACGTACGTA\n")
    assert "|" not in queries.read_text(encoding="utf-8")


def test_collapses_identical_reads_without_abundance_in_query_headers(tmp_path):
    path = tmp_path / "sample.fastq.gz"
    contents = (
        VALID_LIBRARY
        + "@read_three\nACGTACGTACGTACGTACGTA\n+\nIIIIIIIIIIIIIIIIIIIII\n"
    )
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(contents)

    records = read_fastq_collapsed(path)

    assert [
        (record.query_id, record.sequence, record.abundance) for record in records
    ] == [
        ("read_1", "ACGTACGTACGTACGTACGTA", 2),
        ("read_2", "TGCATGCATGCATGCATGCAT", 1),
    ]
    queries = tmp_path / "collapsed_queries.fasta"
    write_mapping_queries(records, queries)
    assert queries.read_text(encoding="utf-8") == (
        ">read_1\nACGTACGTACGTACGTACGTA\n"
        ">read_2\nTGCATGCATGCATGCATGCAT\n"
    )


def test_fastq_ingestion_also_drops_complete_n_records(tmp_path):
    path = tmp_path / "sample.fastq.gz"
    contents = (
        "@ambiguous\nACGTNCGT\n+\nIIIIIIII\n"
        "@retained\nACGTACGT\n+\nIIIIIIII\n"
    )
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(contents)

    records, summary = read_collapsed_fastq(path)

    assert [(record.query_id, record.sequence, record.abundance) for record in records] == [
        ("read_1", "ACGTACGT", 1)
    ]
    assert summary.input_records == 2
    assert summary.excluded_n_records == 1
    assert summary.excluded_n_abundance == 1
    assert summary.retained_abundance == 1


def test_reads_collapsed_tags_drops_complete_n_rows_and_audits_abundance(tmp_path):
    path = tmp_path / "sample.tag"
    path.write_text(
        "ACGT\t5\n"
        "NACGT\t7\n"
        "ACNGT\t11\n"
        "ACGU\t3\n"
        "ACGT\t2\n",
        encoding="utf-8",
    )

    records, summary = read_collapsed_tag(path)

    assert [
        (record.query_id, record.sequence, record.abundance) for record in records
    ] == [("read_1", "ACGT", 10)]
    assert summary.input_records == 5
    assert summary.input_total_abundance == 28
    assert summary.excluded_n_records == 2
    assert summary.excluded_n_abundance == 18
    assert summary.retained_records == 3
    assert summary.retained_unique_sequences == 1
    assert summary.retained_abundance == 10


def test_collapsed_tag_length_filter_is_explicit_and_audited_after_n_exclusion(tmp_path):
    path = tmp_path / "sample.tag"
    path.write_text(
        "ACGT\t5\n"
        "ACGTA\t7\n"
        "ACGTAC\t11\n"
        "ACNGTAC\t13\n",
        encoding="utf-8",
    )

    records, summary = read_collapsed_tag(path, min_length=5, max_length=5)

    assert [(record.sequence, record.abundance) for record in records] == [("ACGTA", 7)]
    assert summary.excluded_n_records == 1
    assert summary.excluded_n_abundance == 13
    assert summary.excluded_length_records == 2
    assert summary.excluded_length_abundance == 16
    assert summary.retained_records == 1
    assert summary.retained_abundance == 7


@pytest.mark.parametrize(
    "contents",
    [
        "ACGT\n",
        "ACGT\tfive\n",
        "ACGT\t0\n",
        "ACRT\t5\n",
    ],
)
def test_rejects_malformed_or_non_n_ambiguous_collapsed_tags(tmp_path, contents):
    path = tmp_path / "invalid.tag"
    path.write_text(contents, encoding="utf-8")

    with pytest.raises(ValidationError):
        read_collapsed_tag(path)


@pytest.mark.parametrize("contents", ["@read\nACGT\n+\n", "@read\nACGT\n-\nIIII\n"])
def test_rejects_malformed_fastq(tmp_path, contents):
    path = tmp_path / "invalid.fastq.gz"
    path.write_bytes(gzip.compress(contents.encode()))
    with pytest.raises(ValidationError):
        read_fastq(path)


def test_only_src_libraries_are_resolved_and_all_libs_is_excluded(tmp_path):
    processed = tmp_path / "processed_libraries"
    processed.mkdir()
    (processed / "Epi_1.fastq.gz").write_bytes(gzip.compress(VALID_LIBRARY.encode()))
    (processed / "ignored.fas.gz").write_bytes(gzip.compress(b">seq_1|1\nACGT\n"))
    targets = tmp_path / "targets.tsv"
    targets.write_text(
        "sample_id\tcondition\tbiological_replicate\tlibrary_id\n"
        "Epi_1\tEpi\t1\tEpi_1\n",
        encoding="utf-8",
    )

    found = find_processed_libraries(processed, library_format="fastq")
    resolved = resolve_target_libraries(
        read_targets(targets), processed, library_format="fastq"
    )

    assert set(found) == {"Epi_1"}
    assert resolved["Epi_1"].name == "Epi_1.fastq.gz"


def test_resolves_collapsed_tag_libraries_by_declared_library_id(tmp_path):
    processed = tmp_path / "processed_libraries"
    processed.mkdir()
    (processed / "WT_0.4_r1.tag").write_text("ACGT\t5\n", encoding="utf-8")
    targets = tmp_path / "targets.tsv"
    targets.write_text(
        "sample_id\tcondition\tbiological_replicate\tlibrary_id\n"
        "WT_0p4_r1\tWT_0p4\t1\tWT_0.4_r1\n",
        encoding="utf-8",
    )

    found = find_processed_libraries(processed, library_format="tag")
    resolved = resolve_target_libraries(
        read_targets(targets), processed, library_format="tag"
    )

    assert set(found) == {"WT_0.4_r1"}
    assert resolved["WT_0p4_r1"].name == "WT_0.4_r1.tag"
