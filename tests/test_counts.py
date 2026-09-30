from __future__ import annotations

from phasis_srna_da.counts import CountAccumulator, Hierarchy
from phasis_srna_da.models import MappingAssignment, ReferenceRecord


def _record(atomic_id, srna_class, feature_id, sequence, **kwargs):
    return ReferenceRecord(
        atomic_id=atomic_id,
        srna_class=srna_class,
        feature_id=feature_id,
        sequence=sequence,
        source_path="fixture",
        source_kind=kwargs.pop("source_kind", "reference_fasta"),
        **kwargs,
    )


def test_equivalence_groups_preserve_abundance_and_phasi_projections_are_conservative():
    tag = "ACGTACGTACGTACGTACGTA"
    records = (
        _record("ref_1", "miRNA", "miR156a", "AAAAAAAAAAAAAAAAAAAAA"),
        _record("ref_2", "miRNA", "miR156b", "AAAAAAAAAAAAAAAAAAAAA"),
        _record("ref_3", "tRNA", "tRNA-Gly-1", "CCCCCCCCCCCCCCCCCCCCC"),
        _record(
            "ref_4",
            "phasiRNA-21",
            "phase21|locusA|strand=w|observed=100|" + tag,
            tag,
            source_kind="phasis_product_catalog",
            phase=21,
            locus_id="locusA",
            tag_sequence=tag,
            strand="w",
            observed_pos="100",
        ),
        _record(
            "ref_5",
            "phasiRNA-21",
            "phase21|locusB|strand=w|observed=200|" + tag,
            tag,
            source_kind="phasis_product_catalog",
            phase=21,
            locus_id="locusB",
            tag_sequence=tag,
            strand="w",
            observed_pos="200",
        ),
    )
    hierarchy = Hierarchy(by_level={"family": {("miRNA", "miR156a"): "miR156", ("miRNA", "miR156b"): "miR156"}})
    accumulator = CountAccumulator(["s1"], records, hierarchy=hierarchy)
    assignments = [
        MappingAssignment("s1", "seq_1", "A" * 21, 50, ("ref_1", "ref_2")),
        MappingAssignment("s1", "seq_2", "C" * 21, 7, ("ref_3",)),
        MappingAssignment("s1", "seq_3", "G" * 21, 8, ()),
        MappingAssignment("s1", "seq_4", tag, 9, ("ref_4", "ref_5")),
        MappingAssignment("s1", "seq_5", tag, 11, ("ref_3", "ref_4")),
    ]

    list(accumulator.add_assignments(assignments))
    primary = accumulator.projections["global_atomic"]
    primary_counts = {unit: counts["s1"] for unit, counts in primary.counts.items()}
    tag_counts = accumulator.projections["phasi_tag_seq"].counts

    assert primary_counts["ref_3"] == 7
    assert sorted(value for unit, value in primary_counts.items() if unit.startswith("EQ__")) == [9, 11, 50]
    assert accumulator.projections["parent_family"].counts["parent::family::miR156"]["s1"] == 50
    assert tag_counts[f"phase21::tag::{tag}"]["s1"] == 9
    assert [
        counts["s1"]
        for unit, counts in accumulator.projections["phasi_tag_seq"].counts.items()
        if unit.startswith("EQ__")
    ] == [11]
    assert len(accumulator.projections["phasi_locus"].equivalence_members) == 2
    conservation = primary.conservation["s1"]
    assert conservation == {
        "input_abundance": 85,
        "unique_assigned_abundance": 7,
        "equivalence_assigned_abundance": 70,
        "unassigned_abundance": 8,
    }
