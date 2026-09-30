"""One-index, exact Bowtie 1 mapping and SAM-to-assignment conversion."""

from __future__ import annotations

import shlex
import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from .errors import ExternalToolError, ValidationError
from .models import FastqRead, MappingAssignment, ReferenceRecord


@dataclass(frozen=True)
class BowtieRun:
    """Command provenance from building an index or mapping a library."""

    command: tuple[str, ...]
    stdout_path: str | None
    stderr_path: str

    @property
    def rendered_command(self) -> str:
        return shlex.join(self.command)


def _resolve_executable(value: str, label: str) -> str:
    candidate = shutil.which(value)
    if candidate is None:
        raise ExternalToolError(
            f"Could not find {label} executable {value!r} on PATH. Install Bowtie 1 "
            "or pass the executable path explicitly. Bowtie 2 is not supported."
        )
    return candidate


def executable_version(value: str) -> str:
    """Best-effort version capture for a command recorded in provenance."""

    executable = _resolve_executable(value, value)
    result = subprocess.run(
        [executable, "--version"], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False
    )
    output = result.stdout.strip()
    return output if output else f"exit={result.returncode}; no version output"


def build_bowtie_index(
    reference_fasta: Path,
    index_prefix: Path,
    *,
    bowtie_build: str = "bowtie-build",
    stderr_path: Path,
) -> BowtieRun:
    """Build a Bowtie 1 index from the atomic competitive reference catalog."""

    executable = _resolve_executable(bowtie_build, "bowtie-build")
    index_prefix.parent.mkdir(parents=True, exist_ok=True)
    command = (executable, str(reference_fasta), str(index_prefix))
    with stderr_path.open("wt", encoding="utf-8") as stderr:
        result = subprocess.run(command, stdout=stderr, stderr=subprocess.STDOUT, check=False)
    if result.returncode != 0:
        raise ExternalToolError(
            f"bowtie-build failed with exit code {result.returncode}; see {stderr_path}."
        )
    return BowtieRun(command=command, stdout_path=None, stderr_path=str(stderr_path))


def run_bowtie_exact(
    index_prefix: Path,
    query_fasta: Path,
    sam_path: Path,
    *,
    threads: int,
    bowtie: str = "bowtie",
    stderr_path: Path,
) -> BowtieRun:
    """Map one canonical library with all direct perfect Bowtie 1 alignments.

    ``-a`` is essential: assignment logic needs every target relationship.
    ``--norc`` prevents a reverse-complement match from being silently treated
    as the same small-RNA sequence. Multiple positions in one target record are
    later deduplicated to one target ID.
    """

    if threads < 1:
        raise ValidationError("--threads must be at least 1.")
    executable = _resolve_executable(bowtie, "bowtie")
    sam_path.parent.mkdir(parents=True, exist_ok=True)
    stderr_path.parent.mkdir(parents=True, exist_ok=True)
    command = (
        executable,
        "-f",
        "-v",
        "0",
        "-a",
        "--best",
        "--strata",
        "--norc",
        "-p",
        str(threads),
        "-S",
        str(index_prefix),
        str(query_fasta),
    )
    with sam_path.open("wt", encoding="utf-8") as stdout, stderr_path.open(
        "wt", encoding="utf-8"
    ) as stderr:
        result = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
    if result.returncode != 0:
        raise ExternalToolError(
            f"bowtie failed with exit code {result.returncode}; see {stderr_path}."
        )
    return BowtieRun(command=command, stdout_path=str(sam_path), stderr_path=str(stderr_path))


def _sam_query_id(raw_name: str) -> str:
    """Validate the internal, abundance-free library query ID in a SAM record."""

    query_id = raw_name
    if not query_id.startswith("read_"):
        raise ValidationError(
            f"SAM QNAME {raw_name!r} is not an internal FASTQ query identifier."
        )
    return query_id


def assignments_from_sam(
    sample_id: str,
    records: list[FastqRead],
    sam_path: Path,
    reference_records: tuple[ReferenceRecord, ...],
) -> list[MappingAssignment]:
    """Parse all alignment relationships while conserving input abundance.

    A query mapping to repeated locations in a single full tRNA record appears
    in several SAM rows but yields one atomic target. A query mapping to several
    atomic references yields a later equivalence group, not several counts.
    Identical reads or tag rows may be represented by one abundance-bearing query;
    abundance is kept outside the query identifier and assigned once.
    """

    if not sam_path.is_file():
        raise ValidationError(f"SAM output does not exist: {sam_path}")
    expected = {record.query_id: record for record in records}
    references_by_id = {record.atomic_id: record for record in reference_records}
    hits: dict[str, set[str]] = defaultdict(set)
    with sam_path.open("rt", encoding="utf-8", newline="") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if raw_line.startswith("@"):
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) < 11:
                raise ValidationError(f"{sam_path}:{line_number}: malformed SAM record.")
            query_id = _sam_query_id(fields[0])
            if query_id not in expected:
                raise ValidationError(
                    f"{sam_path}:{line_number}: QNAME {fields[0]!r} was not present in the input library."
                )
            try:
                flag = int(fields[1])
            except ValueError as exc:
                raise ValidationError(f"{sam_path}:{line_number}: invalid SAM flag {fields[1]!r}.") from exc
            reference_id = fields[2]
            if flag & 4 or reference_id == "*":
                continue
            reference_record = references_by_id.get(reference_id)
            if reference_record is None:
                raise ValidationError(
                    f"{sam_path}:{line_number}: alignment names unknown reference {reference_id!r}."
                )
            # Bowtie allows a 21-nt query to align inside a 24-nt reference.
            # That is appropriate for full structural-RNA records but not for
            # a fixed Phasis product: a 21-nt read is not the same physical
            # product as a 24-nt tag merely because it matches a substring.
            if (
                reference_record.source_kind == "phasis_product_catalog"
                and expected[query_id].sequence != reference_record.sequence
            ):
                continue
            hits[query_id].add(reference_id)

    assignments: list[MappingAssignment] = []
    for record in records:
        assignments.append(
            MappingAssignment(
                sample_id=sample_id,
                query_id=record.query_id,
                sequence=record.sequence,
                abundance=record.abundance,
                atomic_ids=tuple(sorted(hits.get(record.query_id, set()))),
            )
        )
    return assignments
