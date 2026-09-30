"""Build a fixed, non-discovering small-RNA reference catalog for Bowtie 1."""

from __future__ import annotations

import csv
import hashlib
from dataclasses import dataclass
from pathlib import Path

from .errors import ValidationError
from .io import canonical_sequence, parse_reference_spec, read_reference_fasta
from .models import PhasiRNAAssociation, ReferenceRecord


PHASI_REQUIRED_COLUMNS = {
    "identifier",
    "phase",
    "window_unit_id",
    "strand",
    "expected_register_pos",
    "register_class",
    "tag_seq",
}
PHASI_ALLOWED_CLASSES = {"core_exact", "core_offset", "extended_exact"}
PHASI_REFERENCE_PREFIX = "phasiRNA-"
SUPPORTED_PHASES = {21, 24}


@dataclass(frozen=True)
class CatalogDecision:
    """One original Phasis table row and the deterministic inclusion decision."""

    association: PhasiRNAAssociation
    decision: str
    selected: bool


@dataclass(frozen=True)
class Catalog:
    """The atomic mapping reference plus retained Phasis association provenance."""

    records: tuple[ReferenceRecord, ...]
    phasi_decisions: tuple[CatalogDecision, ...]


def _read_phasi_rows(path: Path) -> list[PhasiRNAAssociation]:
    if not path.is_file():
        raise ValidationError(f"Phasis product catalog does not exist: {path}")
    if "phas_like" in path.name.lower():
        raise ValidationError(
            f"{path}: PHAS_like tables are not accepted as a fixed phasiRNA product catalog."
        )
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(PHASI_REQUIRED_COLUMNS - set(reader.fieldnames or []))
        if missing:
            raise ValidationError(
                f"{path}: missing required Phasis columns: {', '.join(missing)}."
            )
        fields_in_source = list(reader.fieldnames or [])
        rows: list[PhasiRNAAssociation] = []
        for row_number, row in enumerate(reader, start=2):
            value = {key: (item or "").strip() for key, item in row.items()}
            missing_values = [key for key in PHASI_REQUIRED_COLUMNS if not value.get(key)]
            if missing_values:
                raise ValidationError(
                    f"{path}:{row_number}: blank required Phasis value(s): "
                    f"{', '.join(sorted(missing_values))}."
                )
            try:
                phase = int(value["phase"])
            except ValueError as exc:
                raise ValidationError(
                    f"{path}:{row_number}: phase must be an integer, got {value['phase']!r}."
                ) from exc
            if phase <= 0:
                raise ValidationError(f"{path}:{row_number}: phase must be positive.")
            if phase not in SUPPORTED_PHASES:
                raise ValidationError(
                    f"{path}:{row_number}: only 21- and 24-nt Phasis product catalogs are "
                    f"supported in v1, got phase={phase}."
                )
            if value["register_class"] not in PHASI_ALLOWED_CLASSES:
                raise ValidationError(
                    f"{path}:{row_number}: unsupported register_class "
                    f"{value['register_class']!r}; expected one of "
                    f"{', '.join(sorted(PHASI_ALLOWED_CLASSES))}."
                )
            tag_sequence = canonical_sequence(
                value["tag_seq"], context=f"{path}:{row_number} tag_seq"
            )
            if len(tag_sequence) != phase:
                raise ValidationError(
                    f"{path}:{row_number}: tag_seq length {len(tag_sequence)} does not match "
                    f"phase {phase}."
                )
            canonical_source_row = "\x1f".join(
                f"{field}={(row.get(field) or '').strip()}" for field in fields_in_source
            )
            rows.append(
                PhasiRNAAssociation(
                    phase=phase,
                    locus_id=value["identifier"],
                    tag_sequence=tag_sequence,
                    register_class=value["register_class"],
                    c_id=value.get("cID", ""),
                    alib=value.get("alib", ""),
                    window_unit_id=value["window_unit_id"],
                    window_unit_role=value.get("window_unit_role", ""),
                    window_unit_rank=value.get("window_unit_rank", ""),
                    window_unit_shift_nt=value.get("window_unit_shift_nt", ""),
                    strand=value["strand"],
                    expected_register_pos=value["expected_register_pos"],
                    observed_pos=value.get("observed_pos", ""),
                    abun=value.get("abun", ""),
                    hits=value.get("hits", ""),
                    source_path=str(path),
                    source_row=row_number,
                    source_row_sha256=hashlib.sha256(
                        canonical_source_row.encode("utf-8")
                    ).hexdigest(),
                )
            )
    if not rows:
        raise ValidationError(f"{path}: no phasiRNA product rows found.")
    return rows


def _core_register_key(row: PhasiRNAAssociation) -> tuple[str | int, ...]:
    """The agreed scope in which a core exact product suppresses an offset call."""

    return (
        row.phase,
        row.locus_id,
        row.strand,
        row.window_unit_id,
        row.expected_register_pos,
    )


def select_phasi_associations(
    paths: list[Path], *, allow_catalog_union: bool = False
) -> list[CatalogDecision]:
    """Apply the agreed core/offset/extension product-catalog policy.

    - every ``core_exact`` row is selected;
    - a ``core_offset`` row is selected only when no ``core_exact`` row exists
      for the same phase, locus, strand, window unit, and expected register;
    - every ``extended_exact`` row is selected.

    Selection never uses ``abun``, ``alib``, ``cID``, or ``hits``. Those are
    run-specific provenance, not counts for the downstream re-count.
    """

    all_rows: list[PhasiRNAAssociation] = []
    catalog_paths_by_phase: dict[int, set[str]] = {}
    for path in paths:
        rows = _read_phasi_rows(path)
        all_rows.extend(rows)
        for phase in {row.phase for row in rows}:
            catalog_paths_by_phase.setdefault(phase, set()).add(str(path))
    repeated_phases = {
        phase: source_paths
        for phase, source_paths in catalog_paths_by_phase.items()
        if len(source_paths) > 1
    }
    if repeated_phases and not allow_catalog_union:
        details = "; ".join(
            f"phase {phase}: {', '.join(sorted(source_paths))}"
            for phase, source_paths in sorted(repeated_phases.items())
        )
        raise ValidationError(
            "Multiple Phasis product catalogs define the same phase. Supply one frozen catalog "
            "per phase, or explicitly pass --allow-catalog-union after documenting the curated "
            f"merge. Conflicts: {details}"
        )
    exact_keys = {
        _core_register_key(row) for row in all_rows if row.register_class == "core_exact"
    }
    decisions: list[CatalogDecision] = []
    for row in all_rows:
        if row.register_class == "core_exact":
            decisions.append(CatalogDecision(row, "selected_core_exact", True))
        elif row.register_class == "core_offset":
            if _core_register_key(row) in exact_keys:
                decisions.append(
                    CatalogDecision(row, "excluded_core_exact_present", False)
                )
            else:
                decisions.append(CatalogDecision(row, "selected_core_offset_fallback", True))
        elif row.register_class == "extended_exact":
            decisions.append(CatalogDecision(row, "selected_extended_exact", True))
        else:
            decisions.append(
                CatalogDecision(row, "excluded_unsupported_register_class", False)
            )
    return decisions


def _new_atomic_id(position: int) -> str:
    return f"ref_{position:09d}"


def build_catalog(
    reference_specs: list[str],
    phasi_catalog_paths: list[Path],
    *,
    allow_catalog_union: bool = False,
) -> Catalog:
    """Construct one atomic, competitive reference catalog.

    User reference FASTAs and Phasis product tags are placed in one index so a
    sequence cannot become falsely unique merely because classes were mapped in
    separate runs. Family/locus hierarchy is deliberately not placed in this
    atomic index; it is a later projection of the mapping relationships.
    """

    if not reference_specs:
        raise ValidationError("At least one --reference CLASS=FASTA declaration is required.")
    if not phasi_catalog_paths:
        raise ValidationError("At least one --phasiRNAs *_phasiRNAs.tsv product file is required.")

    basic_records: list[tuple[str, str, str, str]] = []
    seen_feature_keys: set[tuple[str, str]] = set()
    for raw_spec in reference_specs:
        srna_class, path = parse_reference_spec(raw_spec)
        if srna_class.startswith(PHASI_REFERENCE_PREFIX):
            raise ValidationError(
                f"Reference class {srna_class!r} is reserved for products built from --phasiRNAs."
            )
        for feature_id, sequence in read_reference_fasta(path):
            key = (srna_class, feature_id)
            if key in seen_feature_keys:
                raise ValidationError(
                    f"Duplicate reference feature {srna_class}::{feature_id}; each class/ID pair "
                    "must be unique across all --reference files."
                )
            seen_feature_keys.add(key)
            basic_records.append((srna_class, feature_id, sequence, str(path)))

    decisions = select_phasi_associations(
        phasi_catalog_paths, allow_catalog_union=allow_catalog_union
    )
    product_keys: dict[tuple[int, str, str, str, str], list[CatalogDecision]] = {}
    for decision in decisions:
        if not decision.selected:
            continue
        association = decision.association
        # A physical product is strand/observed position/sequence at a Phasis
        # locus. The same product can receive different register annotations in
        # different Phasis runs; those associations must not create copies.
        # When an older table omits observed_pos, expected position is the only
        # available fallback and is retained transparently in the feature ID.
        observed_position = association.observed_pos or association.expected_register_pos
        key = (
            association.phase,
            association.locus_id,
            association.strand,
            observed_position,
            association.tag_sequence,
        )
        product_keys.setdefault(key, []).append(decision)

    if not product_keys:
        raise ValidationError(
            "No countable Phasis products remain after applying the core/offset/extension policy."
        )

    records: list[ReferenceRecord] = []
    for srna_class, feature_id, sequence, source_path in sorted(basic_records):
        records.append(
            ReferenceRecord(
                atomic_id=_new_atomic_id(len(records) + 1),
                srna_class=srna_class,
                feature_id=feature_id,
                sequence=sequence,
                source_path=source_path,
                source_kind="reference_fasta",
            )
        )
    for (
        phase,
        locus_id,
        strand,
        observed_position,
        tag_sequence,
    ), product_decisions in sorted(product_keys.items()):
        source_paths = sorted({item.association.source_path for item in product_decisions})
        feature_id = (
            f"phase{phase}|{locus_id}|strand={strand}|observed={observed_position}|{tag_sequence}"
        )
        records.append(
            ReferenceRecord(
                atomic_id=_new_atomic_id(len(records) + 1),
                srna_class=f"{PHASI_REFERENCE_PREFIX}{phase}",
                feature_id=feature_id,
                sequence=tag_sequence,
                source_path=";".join(source_paths),
                source_kind="phasis_product_catalog",
                phase=phase,
                locus_id=locus_id,
                tag_sequence=tag_sequence,
                strand=strand,
                observed_pos=observed_position,
            )
        )
    if not records:
        raise ValidationError("The supplied references and selected Phasis product catalog are empty.")
    return Catalog(records=tuple(records), phasi_decisions=tuple(decisions))


def write_atomic_reference(records: tuple[ReferenceRecord, ...], destination: Path) -> None:
    """Write the reference that Bowtie indexes, using safe internal IDs only."""

    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("wt", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(f">{record.atomic_id}\n{record.sequence}\n")


def catalog_record_rows(records: tuple[ReferenceRecord, ...]) -> list[dict[str, object]]:
    """Return the atomic-reference manifest rows for a TSV output."""

    return [
        {
            "atomic_id": record.atomic_id,
            "srna_class": record.srna_class,
            "feature_id": record.feature_id,
            "sequence": record.sequence,
            "source_path": record.source_path,
            "source_kind": record.source_kind,
            "phase": record.phase or "",
            "locus_id": record.locus_id or "",
            "tag_sequence": record.tag_sequence or "",
            "strand": record.strand or "",
            "observed_pos": record.observed_pos or "",
        }
        for record in records
    ]


def catalog_decision_rows(catalog: Catalog) -> list[dict[str, object]]:
    """Return every catalog row decision with physical-product provenance.

    Window/register fields remain a many-to-one association ledger. The atomic
    reference ID makes it explicit when several Phasis rows describe one
    physical product rather than independent countable features.
    """

    atomic_by_product_key = {
        (
            record.phase,
            record.locus_id,
            record.strand,
            record.observed_pos,
            record.tag_sequence,
        ): record.atomic_id
        for record in catalog.records
        if record.source_kind == "phasis_product_catalog"
    }
    rows: list[dict[str, object]] = []
    for item in catalog.phasi_decisions:
        association = item.association
        observed_position = association.observed_pos or association.expected_register_pos
        product_key = (
            association.phase,
            association.locus_id,
            association.strand,
            observed_position,
            association.tag_sequence,
        )
        rows.append(
            {
                "selected": str(item.selected).lower(),
                "decision": item.decision,
                "phase": association.phase,
                "identifier": association.locus_id,
                "tag_seq": association.tag_sequence,
                "register_class": association.register_class,
                "cID": association.c_id,
                "alib": association.alib,
                "window_unit_id": association.window_unit_id,
                "window_unit_role": association.window_unit_role,
                "window_unit_rank": association.window_unit_rank,
                "window_unit_shift_nt": association.window_unit_shift_nt,
                "strand": association.strand,
                "expected_register_pos": association.expected_register_pos,
                "observed_pos": association.observed_pos,
                "abun": association.abun,
                "hits": association.hits,
                "source_path": association.source_path,
                "source_row": association.source_row,
                "source_row_sha256": association.source_row_sha256,
                "physical_product_key": "|".join(str(part) for part in product_key),
                "atomic_id": atomic_by_product_key.get(product_key, ""),
            }
        )
    return rows
