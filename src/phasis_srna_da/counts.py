"""Count-conserving projections of exact Bowtie relationships into DA matrices."""

from __future__ import annotations

import csv
import gzip
import hashlib
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator

from .errors import ValidationError
from .models import MappingAssignment, ReferenceRecord


@dataclass(frozen=True)
class Hierarchy:
    """Optional declared parent relationships; names are never inferred."""

    by_level: dict[str, dict[tuple[str, str], str]]


def read_hierarchy(path: Path) -> Hierarchy:
    """Read optional ``class, feature_id, parent_level, parent_id`` TSV metadata."""

    if not path.is_file():
        raise ValidationError(f"Hierarchy catalog does not exist: {path}")
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"class", "feature_id", "parent_level", "parent_id"}
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValidationError(
                f"{path}: hierarchy catalog is missing columns: {', '.join(missing)}."
            )
        by_level: dict[str, dict[tuple[str, str], str]] = defaultdict(dict)
        for row_number, row in enumerate(reader, start=2):
            values = {key: (value or "").strip() for key, value in row.items()}
            missing_values = [key for key in required if not values[key]]
            if missing_values:
                raise ValidationError(
                    f"{path}:{row_number}: blank hierarchy field(s): {', '.join(sorted(missing_values))}."
                )
            key = (values["class"], values["feature_id"])
            level = values["parent_level"]
            parent_id = values["parent_id"]
            existing = by_level[level].get(key)
            if existing is not None and existing != parent_id:
                raise ValidationError(
                    f"{path}:{row_number}: conflicting {level!r} parent for "
                    f"{key[0]}::{key[1]} ({existing!r} versus {parent_id!r})."
                )
            by_level[level][key] = parent_id
    if not by_level:
        raise ValidationError(f"{path}: hierarchy catalog has no rows.")
    return Hierarchy(by_level=dict(by_level))


@dataclass
class ProjectionAccumulator:
    """A raw integer count matrix plus conservation and equivalence evidence."""

    name: str
    sample_ids: tuple[str, ...]
    counts: dict[str, dict[str, int]] = field(default_factory=lambda: defaultdict(lambda: defaultdict(int)))
    unit_metadata: dict[str, dict[str, str]] = field(default_factory=dict)
    equivalence_members: dict[str, tuple[tuple[str, str], ...]] = field(default_factory=dict)
    conservation: dict[str, dict[str, int]] = field(
        default_factory=lambda: defaultdict(
            lambda: {
                "input_abundance": 0,
                "unique_assigned_abundance": 0,
                "equivalence_assigned_abundance": 0,
                "unassigned_abundance": 0,
            }
        )
    )

    def _equivalence_id(self, targets: tuple[tuple[str, str], ...]) -> str:
        # IDs emitted by Bowtie are internal and depend on reference ordering.
        # Hash the canonical biological labels instead so an unchanged mapping
        # relationship gets the same equivalence ID after unrelated references
        # are added to the catalog.
        joined = "\x1f".join(sorted({label for _, label in targets}))
        digest = hashlib.sha256(joined.encode("utf-8")).hexdigest()[:20]
        group_id = f"EQ__{self.name}__{digest}"
        existing = self.equivalence_members.get(group_id)
        if existing is not None and existing != targets:
            raise RuntimeError("Equivalence-group hash collision; please report this bug.")
        self.equivalence_members[group_id] = targets
        return group_id

    def add(
        self,
        *,
        sample_id: str,
        abundance: int,
        targets: Iterable[tuple[str, str]],
        unique_metadata: dict[str, str] | None = None,
    ) -> tuple[str, str]:
        """Assign abundance once to one target, one EQ group, or unassigned."""

        if sample_id not in self.sample_ids:
            raise ValidationError(f"Unknown sample {sample_id!r} in projection {self.name!r}.")
        if abundance <= 0:
            raise ValidationError("Read abundance must be positive.")
        target_tuple = tuple(sorted(set(targets)))
        summary = self.conservation[sample_id]
        summary["input_abundance"] += abundance
        if not target_tuple:
            summary["unassigned_abundance"] += abundance
            return "unassigned", ""
        if len(target_tuple) == 1:
            unit_id, label = target_tuple[0]
            self.counts[unit_id][sample_id] += abundance
            metadata = {
                "unit_id": unit_id,
                "unit_label": label,
                "assignment_type": "unique",
                "projection": self.name,
            }
            if unique_metadata:
                metadata.update(unique_metadata)
            existing = self.unit_metadata.get(unit_id)
            if existing is not None and existing != metadata:
                raise RuntimeError(f"Conflicting metadata for count unit {unit_id!r}.")
            self.unit_metadata[unit_id] = metadata
            summary["unique_assigned_abundance"] += abundance
            return "unique", unit_id
        group_id = self._equivalence_id(target_tuple)
        self.counts[group_id][sample_id] += abundance
        self.unit_metadata.setdefault(
            group_id,
            {
                "unit_id": group_id,
                "unit_label": " | ".join(label for _, label in target_tuple),
                "assignment_type": "equivalence_group",
                "projection": self.name,
            },
        )
        summary["equivalence_assigned_abundance"] += abundance
        return "equivalence_group", group_id


class CountAccumulator:
    """Conservative primary and Phasis-derived count projections.

    Every projection independently counts an input abundance once. Derived
    Phasis projections are intentionally *not* additive to the primary matrix.
    They answer different questions: tag sequence, physical product, locus/tag,
    and locus output.
    """

    def __init__(
        self,
        sample_ids: Iterable[str],
        reference_records: tuple[ReferenceRecord, ...],
        hierarchy: Hierarchy | None = None,
    ) -> None:
        self.sample_ids = tuple(sample_ids)
        self.records = {record.atomic_id: record for record in reference_records}
        if len(self.records) != len(reference_records):
            raise ValidationError("Atomic reference IDs must be unique.")
        self.hierarchy = hierarchy
        self.projections: dict[str, ProjectionAccumulator] = {
            "global_atomic": ProjectionAccumulator("global_atomic", self.sample_ids),
            "phasi_tag_seq": ProjectionAccumulator("phasi_tag_seq", self.sample_ids),
            "phasi_product": ProjectionAccumulator("phasi_product", self.sample_ids),
            "phasi_locus_tag": ProjectionAccumulator("phasi_locus_tag", self.sample_ids),
            "phasi_locus": ProjectionAccumulator("phasi_locus", self.sample_ids),
        }
        if hierarchy is not None:
            for level in sorted(hierarchy.by_level):
                name = f"parent_{level}"
                self.projections[name] = ProjectionAccumulator(name, self.sample_ids)

    @staticmethod
    def _feature_label(record: ReferenceRecord) -> str:
        return f"{record.srna_class}::{record.feature_id}"

    @staticmethod
    def _is_phasi(record: ReferenceRecord) -> bool:
        return record.source_kind == "phasis_product_catalog"

    def _records_for(self, assignment: MappingAssignment) -> list[ReferenceRecord]:
        try:
            return [self.records[atomic_id] for atomic_id in assignment.atomic_ids]
        except KeyError as exc:
            raise ValidationError(
                f"Mapping assignment contains unknown atomic reference {exc.args[0]!r}."
            ) from exc

    def _add_primary(self, assignment: MappingAssignment, records: list[ReferenceRecord]) -> tuple[str, str]:
        targets = [(record.atomic_id, self._feature_label(record)) for record in records]
        metadata = None
        if len(records) == 1:
            record = records[0]
            metadata = {
                "srna_class": record.srna_class,
                "feature_id": record.feature_id,
                "source_kind": record.source_kind,
            }
        return self.projections["global_atomic"].add(
            sample_id=assignment.sample_id,
            abundance=assignment.abundance,
            targets=targets,
            unique_metadata=metadata,
        )

    def _add_phasi_projections(
        self, assignment: MappingAssignment, records: list[ReferenceRecord]
    ) -> None:
        phasi_records = [record for record in records if self._is_phasi(record)]
        non_phasi_records = [record for record in records if not self._is_phasi(record)]
        if not phasi_records:
            for name in ("phasi_tag_seq", "phasi_product", "phasi_locus_tag", "phasi_locus"):
                self.projections[name].add(
                    sample_id=assignment.sample_id,
                    abundance=assignment.abundance,
                    targets=(),
                )
            return
        phases = {record.phase for record in phasi_records}
        tags = {record.tag_sequence for record in phasi_records}
        if len(tags) != 1:
            raise RuntimeError("Exact mappings of one query to different canonical tag sequences are impossible.")
        tag = next(iter(tags))
        cross_class_targets = [
            (f"atomic::{record.atomic_id}", self._feature_label(record))
            for record in non_phasi_records
        ]
        tag_targets = [
            (f"phase{phase}::tag::{tag}", f"phase {phase} tag_seq {tag}")
            for phase in sorted(phases)
        ] + cross_class_targets
        self.projections["phasi_tag_seq"].add(
            sample_id=assignment.sample_id,
            abundance=assignment.abundance,
            targets=tag_targets,
            unique_metadata=(
                {
                    "phase": str(next(iter(phases))),
                    "tag_seq": tag,
                    "srna_class": f"phasiRNA-{next(iter(phases))}",
                }
                if len(phases) == 1 and not non_phasi_records
                else None
            ),
        )

        product_targets = [
            (record.atomic_id, self._feature_label(record)) for record in phasi_records
        ] + cross_class_targets
        product_metadata = None
        if len(phasi_records) == 1 and not non_phasi_records:
            record = phasi_records[0]
            product_metadata = {
                "phase": str(record.phase),
                "identifier": record.locus_id or "",
                "tag_seq": record.tag_sequence or "",
                "strand": record.strand or "",
                "observed_pos": record.observed_pos or "",
                "srna_class": record.srna_class,
            }
        self.projections["phasi_product"].add(
            sample_id=assignment.sample_id,
            abundance=assignment.abundance,
            targets=product_targets,
            unique_metadata=product_metadata,
        )

        locus_tag_targets: dict[tuple[int | None, str | None, str | None], tuple[str, str]] = {}
        locus_targets: dict[tuple[int | None, str | None], tuple[str, str]] = {}
        for record in phasi_records:
            locus_tag_key = (record.phase, record.locus_id, record.tag_sequence)
            locus_tag_targets[locus_tag_key] = (
                f"phase{record.phase}::locus_tag::{record.locus_id}::{record.tag_sequence}",
                f"phase {record.phase} locus {record.locus_id} tag_seq {record.tag_sequence}",
            )
            locus_key = (record.phase, record.locus_id)
            locus_targets[locus_key] = (
                f"phase{record.phase}::locus::{record.locus_id}",
                f"phase {record.phase} locus {record.locus_id}",
            )
        locus_tag_projection_targets = list(locus_tag_targets.values()) + cross_class_targets
        locus_projection_targets = list(locus_targets.values()) + cross_class_targets
        self.projections["phasi_locus_tag"].add(
            sample_id=assignment.sample_id,
            abundance=assignment.abundance,
            targets=locus_tag_projection_targets,
            unique_metadata=(
                {
                    "phase": str(phasi_records[0].phase),
                    "identifier": phasi_records[0].locus_id or "",
                    "tag_seq": tag,
                    "srna_class": f"phasiRNA-{phasi_records[0].phase}",
                }
                if len(locus_tag_targets) == 1 and not non_phasi_records
                else None
            ),
        )
        self.projections["phasi_locus"].add(
            sample_id=assignment.sample_id,
            abundance=assignment.abundance,
            targets=locus_projection_targets,
            unique_metadata=(
                {
                    "phase": str(phasi_records[0].phase),
                    "identifier": phasi_records[0].locus_id or "",
                    "srna_class": f"phasiRNA-{phasi_records[0].phase}",
                }
                if len(locus_targets) == 1 and not non_phasi_records
                else None
            ),
        )

    def _add_parent_projections(
        self, assignment: MappingAssignment, records: list[ReferenceRecord]
    ) -> None:
        if self.hierarchy is None:
            return
        for level, parent_map in self.hierarchy.by_level.items():
            projection = self.projections[f"parent_{level}"]
            targets: list[tuple[str, str]] = []
            all_have_parent = bool(records)
            for record in records:
                parent_id = parent_map.get((record.srna_class, record.feature_id))
                if parent_id is None:
                    all_have_parent = False
                    break
                targets.append((f"parent::{level}::{parent_id}", f"{level} {parent_id}"))
            if not all_have_parent:
                # No declared parent relation means no parent-level claim. The
                # primary global matrix still preserves the exact EQ evidence.
                projection.add(
                    sample_id=assignment.sample_id,
                    abundance=assignment.abundance,
                    targets=(),
                )
                continue
            projection.add(
                sample_id=assignment.sample_id,
                abundance=assignment.abundance,
                targets=targets,
                unique_metadata=(
                    {"parent_level": level, "parent_id": targets[0][0].rsplit("::", 1)[-1]}
                    if len(set(targets)) == 1
                    else None
                ),
            )

    def add_assignments(self, assignments: Iterable[MappingAssignment]) -> Iterator[dict[str, str]]:
        """Update all projections and yield primary-ledger rows for writing."""

        for assignment in assignments:
            records = self._records_for(assignment)
            status, unit_id = self._add_primary(assignment, records)
            self._add_phasi_projections(assignment, records)
            self._add_parent_projections(assignment, records)
            yield {
                "sample_id": assignment.sample_id,
                "query_id": assignment.query_id,
                "sequence": assignment.sequence,
                "abundance": str(assignment.abundance),
                "atomic_ids": ";".join(assignment.atomic_ids),
                "atomic_features": ";".join(self._feature_label(record) for record in records),
                "primary_assignment_type": status,
                "primary_unit_id": unit_id,
            }

    def count_matrix_rows(self, projection_name: str) -> tuple[list[str], list[dict[str, object]]]:
        """Build stable matrix rows for one named projection."""

        projection = self.projections[projection_name]
        metadata_columns = [
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
        ]
        rows: list[dict[str, object]] = []
        for unit_id in sorted(projection.counts):
            row: dict[str, object] = {name: "" for name in metadata_columns}
            row.update(projection.unit_metadata[unit_id])
            row.update({sample_id: projection.counts[unit_id].get(sample_id, 0) for sample_id in self.sample_ids})
            rows.append(row)
        return metadata_columns + list(self.sample_ids), rows

    def equivalence_rows(self) -> tuple[list[str], list[dict[str, object]]]:
        """Return all equivalence membership rows for every projection."""

        rows: list[dict[str, object]] = []
        for projection_name, projection in sorted(self.projections.items()):
            for group_id, members in sorted(projection.equivalence_members.items()):
                for member_index, (target_id, target_label) in enumerate(members, start=1):
                    rows.append(
                        {
                            "projection": projection_name,
                            "equivalence_group_id": group_id,
                            "member_rank": member_index,
                            "target_id": target_id,
                            "target_label": target_label,
                        }
                    )
        return (
            ["projection", "equivalence_group_id", "member_rank", "target_id", "target_label"],
            rows,
        )

    def conservation_rows(self) -> tuple[list[str], list[dict[str, object]]]:
        """Return per-sample count-conservation summaries for all projections."""

        rows: list[dict[str, object]] = []
        for projection_name, projection in sorted(self.projections.items()):
            for sample_id in self.sample_ids:
                values = projection.conservation[sample_id]
                assigned = values["unique_assigned_abundance"] + values["equivalence_assigned_abundance"]
                reconciled = assigned + values["unassigned_abundance"]
                if values["input_abundance"] != reconciled:
                    raise RuntimeError(
                        f"Count conservation failed for {projection_name}/{sample_id}: "
                        f"{values['input_abundance']} != {reconciled}."
                    )
                rows.append(
                    {
                        "projection": projection_name,
                        "sample_id": sample_id,
                        **values,
                        "assigned_abundance": assigned,
                        "conserved": "true",
                    }
                )
        return (
            [
                "projection",
                "sample_id",
                "input_abundance",
                "unique_assigned_abundance",
                "equivalence_assigned_abundance",
                "unassigned_abundance",
                "assigned_abundance",
                "conserved",
            ],
            rows,
        )


def write_gzip_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, object]]) -> None:
    """Write a gzip-compressed audit table without retaining it in memory."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})
