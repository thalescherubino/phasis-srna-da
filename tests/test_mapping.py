from __future__ import annotations

from phasis_srna_da.mapping import assignments_from_sam
from phasis_srna_da.models import FastqRead, ReferenceRecord


def _reference(atomic_id):
    return ReferenceRecord(
        atomic_id=atomic_id,
        srna_class="tRNA",
        feature_id="tRNA-Gly-1",
        sequence="A" * 21,
        source_path="fixture",
        source_kind="reference_fasta",
    )


def _phasi_reference(atomic_id, sequence, phase):
    return ReferenceRecord(
        atomic_id=atomic_id,
        srna_class=f"phasiRNA-{phase}",
        feature_id=f"phase{phase}-product",
        sequence=sequence,
        source_path="fixture",
        source_kind="phasis_product_catalog",
        phase=phase,
        locus_id="1:100..200",
        tag_sequence=sequence,
        strand="w",
        observed_pos="100",
    )


def test_repeated_alignments_to_one_reference_are_deduplicated(tmp_path):
    sam = tmp_path / "library.sam"
    sam.write_text(
        "@HD\tVN:1.0\n"
        "read_1\t0\tref_1\t1\t255\t21M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAA\t*\n"
        "read_1\t0\tref_1\t20\t255\t21M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAA\t*\n"
        "read_2\t4\t*\t0\t0\t*\t*\t0\t0\tCCCCCCCCCCCCCCCCCCCCC\t*\n",
        encoding="utf-8",
    )
    records = [
        FastqRead("read_1", "A" * 21, abundance=7),
        FastqRead("read_2", "C" * 21),
    ]

    assignments = assignments_from_sam("s1", records, sam, (_reference("ref_1"),))

    assert assignments[0].atomic_ids == ("ref_1",)
    assert assignments[0].abundance == 7
    assert assignments[1].atomic_ids == ()


def test_short_read_is_not_assigned_to_a_longer_phasi_product_substring(tmp_path):
    query = "ACGTACGTACGTACGTACGTA"
    longer_product = query + "CGT"
    sam = tmp_path / "library.sam"
    sam.write_text(
        "@HD\tVN:1.0\n"
        f"read_1\t0\tref_21\t1\t255\t21M\t*\t0\t0\t{query}\t*\n"
        f"read_1\t0\tref_24\t1\t255\t21M\t*\t0\t0\t{query}\t*\n",
        encoding="utf-8",
    )

    assignments = assignments_from_sam(
        "s1",
        [FastqRead("read_1", query)],
        sam,
        (_phasi_reference("ref_21", query, 21), _phasi_reference("ref_24", longer_product, 24)),
    )

    assert assignments[0].atomic_ids == ("ref_21",)
