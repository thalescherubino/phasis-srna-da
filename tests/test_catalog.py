from __future__ import annotations

from phasis_srna_da.catalog import build_catalog, catalog_decision_rows


HEADER = (
    "identifier\tcID\talib\tphase\twindow_unit_id\twindow_unit_role\twindow_unit_rank\t"
    "window_unit_shift_nt\tstrand\tobserved_pos\texpected_register_pos\tregister_class\tabun\ttag_seq\thits\n"
)


def _row(tag, observed, expected, register_class, alib="s1", unit="unit_main"):
    return "\t".join(
        [
            "1:100..200",
            "unused-cid",
            alib,
            "21",
            unit,
            "main",
            "1",
            "0",
            "w",
            str(observed),
            str(expected),
            register_class,
            "999",
            tag,
            "1",
        ]
    ) + "\n"


def test_catalog_uses_exact_offset_fallback_extension_and_physical_product_dedup(tmp_path):
    reference = tmp_path / "mirna.fa"
    reference.write_text(">miR156a\nAAAAAAAAAAAAAAAAAAAAA\n", encoding="utf-8")
    catalog = tmp_path / "21_phasiRNAs.tsv"
    tag_a = "AAAAAAAAAAAAAAAAAAAAA"
    tag_b = "CCCCCCCCCCCCCCCCCCCCC"
    tag_c = "GGGGGGGGGGGGGGGGGGGGG"
    tag_d = "TTTTTTTTTTTTTTTTTTTTT"
    tag_e = "ACGTACGTACGTACGTACGTA"
    tag_f = "CGTACGTACGTACGTACGTAC"
    catalog.write_text(
        HEADER
        + _row(tag_a, 121, 121, "core_exact", "s1")
        + _row(tag_a, 121, 121, "extended_exact", "s2")
        + _row(tag_b, 141, 142, "core_offset", "s1")
        + _row(tag_c, 142, 142, "core_exact", "s2")
        + _row(tag_d, 162, 163, "core_offset", "s1")
        + _row(tag_e, 184, 184, "extended_exact", "s1", "unit_secondary_1")
        + _row(tag_f, 205, 205, "core_exact", "s1")
        + _row(tag_f, 205, 204, "core_offset", "s1"),
        encoding="utf-8",
    )

    result = build_catalog([f"miRNA={reference}"], [catalog])
    decisions = {item.association.tag_sequence + item.association.register_class + item.association.expected_register_pos: item for item in result.phasi_decisions}
    products = [record for record in result.records if record.source_kind == "phasis_product_catalog"]

    assert decisions[tag_b + "core_offset" + "142"].selected is False
    assert decisions[tag_b + "core_offset" + "142"].decision == "excluded_core_exact_present"
    assert decisions[tag_d + "core_offset" + "163"].selected is True
    assert decisions[tag_e + "extended_exact" + "184"].selected is True
    # tag_a has exact+extension memberships, and tag_f has exact+offset memberships;
    # each is one physical reference product rather than multiple counted copies.
    assert len(products) == 5
    assert len([record for record in products if record.tag_sequence == tag_a]) == 1
    assert len([record for record in products if record.tag_sequence == tag_f]) == 1
    association_rows = catalog_decision_rows(result)
    tag_a_rows = [row for row in association_rows if row["tag_seq"] == tag_a]
    assert {row["atomic_id"] for row in tag_a_rows} == {
        next(record.atomic_id for record in products if record.tag_sequence == tag_a)
    }
    assert all(row["source_row_sha256"] for row in association_rows)
