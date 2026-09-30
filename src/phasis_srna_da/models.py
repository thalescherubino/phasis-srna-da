"""Typed records shared by parsing, reference construction, and counting."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Sample:
    """One biological library declared in the user-supplied target file."""

    sample_id: str
    condition: str
    biological_replicate: str
    library_id: str
    batch: str | None = None
    technical_replicate: str | None = None


@dataclass(frozen=True)
class FastqRead:
    """One retained sequence query with abundance kept outside its identifier."""

    query_id: str
    sequence: str
    abundance: int = 1


@dataclass(frozen=True)
class ReferenceRecord:
    """An atomic reference record used in the one competitive Bowtie index."""

    atomic_id: str
    srna_class: str
    feature_id: str
    sequence: str
    source_path: str
    source_kind: str
    phase: int | None = None
    locus_id: str | None = None
    tag_sequence: str | None = None
    strand: str | None = None
    observed_pos: str | None = None


@dataclass(frozen=True)
class PhasiRNAAssociation:
    """Provenance retained for a selected Phasis product catalog row."""

    phase: int
    locus_id: str
    tag_sequence: str
    register_class: str
    c_id: str
    alib: str
    window_unit_id: str
    window_unit_role: str
    window_unit_rank: str
    window_unit_shift_nt: str
    strand: str
    expected_register_pos: str
    observed_pos: str
    abun: str
    hits: str
    source_path: str
    source_row: int
    source_row_sha256: str


@dataclass(frozen=True)
class MappingAssignment:
    """All atomic reference records matched by one retained sequence query."""

    sample_id: str
    query_id: str
    sequence: str
    abundance: int
    atomic_ids: tuple[str, ...]
