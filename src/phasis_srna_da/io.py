"""Strict readers and writers for quality-controlled FASTQ and simple TSV/FASTA inputs."""

from __future__ import annotations

import csv
import gzip
import hashlib
import re
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from .errors import ValidationError
from .models import FastqRead, Sample


SAFE_SEQUENCE = re.compile(r"^[ACGT]+$")
SAFE_CLASS = re.compile(r"^[A-Za-z][A-Za-z0-9_-]*$")
LIBRARY_SUFFIXES = {"fastq": ".fastq.gz", "tag": ".tag"}


@dataclass(frozen=True)
class LibraryReadSummary:
    """Audited ingestion totals before and after whole-record N exclusion."""

    input_records: int
    input_total_abundance: int
    excluded_n_records: int
    excluded_n_abundance: int
    excluded_length_records: int
    excluded_length_abundance: int
    retained_records: int
    retained_unique_sequences: int
    retained_abundance: int


def canonical_sequence(value: str, *, context: str) -> str:
    """Return an uppercase DNA-form sequence, treating RNA U as T.

    The conversion is explicit and is applied equally to query and reference
    sequences before an exact Bowtie search. Ambiguous bases are rejected:
    an exact assignment should not depend on aligner-specific N handling.
    """

    sequence = "".join(value.split()).upper().replace("U", "T")
    if not sequence or not SAFE_SEQUENCE.fullmatch(sequence):
        raise ValidationError(
            f"{context}: expected a non-empty A/C/G/T/U sequence; got {value!r}."
        )
    return sequence


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("rt", encoding="utf-8", newline="")


def _normalized_input_sequence(value: str) -> str:
    """Normalize case and whitespace before applying the explicit N policy."""

    return "".join(value.split()).upper()


def _validate_length_bounds(
    min_length: int | None, max_length: int | None
) -> None:
    if min_length is not None and min_length < 1:
        raise ValidationError("Minimum library sequence length must be at least 1.")
    if max_length is not None and max_length < 1:
        raise ValidationError("Maximum library sequence length must be at least 1.")
    if (
        min_length is not None
        and max_length is not None
        and min_length > max_length
    ):
        raise ValidationError("Minimum library sequence length exceeds maximum length.")


def _outside_length_bounds(
    sequence: str, min_length: int | None, max_length: int | None
) -> bool:
    return (
        (min_length is not None and len(sequence) < min_length)
        or (max_length is not None and len(sequence) > max_length)
    )


def _iter_fastq_rows(path: Path) -> Iterator[tuple[int, str]]:
    """Yield validated FASTQ record numbers and uncanonicalized sequences."""

    if not path.name.endswith(".fastq.gz"):
        raise ValidationError(f"{path}: expected a quality-controlled .fastq.gz file.")
    if not path.is_file():
        raise ValidationError(f"FASTQ library does not exist: {path}")

    with _open_text(path) as handle:
        record_number = 0
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            separator = handle.readline()
            quality = handle.readline()
            record_number += 1
            if not sequence or not separator or not quality:
                raise ValidationError(f"{path}: truncated FASTQ record {record_number}.")
            header = header.rstrip("\r\n")
            sequence = sequence.rstrip("\r\n")
            separator = separator.rstrip("\r\n")
            quality = quality.rstrip("\r\n")
            if not header.startswith("@"):
                raise ValidationError(
                    f"{path}: FASTQ record {record_number} header does not start with '@'."
                )
            if not separator.startswith("+"):
                raise ValidationError(
                    f"{path}: FASTQ record {record_number} separator does not start with '+'."
                )
            if len(sequence) != len(quality):
                raise ValidationError(
                    f"{path}: FASTQ record {record_number} sequence and quality lengths differ."
                )
            yield record_number, sequence
    if record_number == 0:
        raise ValidationError(f"{path}: no FASTQ records found.")


def iter_fastq(path: Path) -> Iterator[FastqRead]:
    """Parse a quality-controlled four-line ``.fastq.gz`` library.

    The original FASTQ headers are not used as Bowtie query IDs. Stable internal
    ``read_N`` IDs make all retained reads count once without adding abundance
    information to a header.
    """

    retained_number = 0
    for record_number, raw_sequence in _iter_fastq_rows(path):
        normalized = _normalized_input_sequence(raw_sequence)
        if "N" in normalized:
            continue
        retained_number += 1
        yield FastqRead(
            query_id=f"read_{retained_number}",
            sequence=canonical_sequence(
                normalized, context=f"{path}: FASTQ record {record_number}"
            ),
        )
    if retained_number == 0:
        raise ValidationError(f"{path}: no records remain after dropping sequences containing N.")


def read_fastq(path: Path) -> list[FastqRead]:
    """Materialize a quality-controlled FASTQ library after full validation."""

    return list(iter_fastq(path))


def read_fastq_collapsed(path: Path) -> list[FastqRead]:
    """Validate and collapse identical reads without abundance in query IDs.

    The mapping FASTA still uses stable, abundance-free read_N identifiers.
    Multiplicity is carried separately on FastqRead and is restored during
    assignment, so every original FASTQ record contributes one count.
    """

    records, _ = read_collapsed_fastq(path)
    return records


def _records_from_abundances(abundance_by_sequence: dict[str, int]) -> list[FastqRead]:
    """Create stable, abundance-free query IDs from first-seen sequence order."""

    return [
        FastqRead(
            query_id=f"read_{index}",
            sequence=sequence,
            abundance=abundance,
        )
        for index, (sequence, abundance) in enumerate(
            abundance_by_sequence.items(), start=1
        )
    ]


def read_collapsed_fastq(
    path: Path,
    *,
    min_length: int | None = None,
    max_length: int | None = None,
) -> tuple[list[FastqRead], LibraryReadSummary]:
    """Collapse FASTQ reads and audit N and explicit length exclusions."""

    _validate_length_bounds(min_length, max_length)
    abundance_by_sequence: dict[str, int] = {}
    input_records = 0
    excluded_n_records = 0
    excluded_length_records = 0
    for record_number, raw_sequence in _iter_fastq_rows(path):
        input_records += 1
        normalized = _normalized_input_sequence(raw_sequence)
        if "N" in normalized:
            excluded_n_records += 1
            continue
        sequence = canonical_sequence(
            normalized, context=f"{path}: FASTQ record {record_number}"
        )
        if _outside_length_bounds(sequence, min_length, max_length):
            excluded_length_records += 1
            continue
        abundance_by_sequence[sequence] = abundance_by_sequence.get(sequence, 0) + 1
    retained_records = input_records - excluded_n_records - excluded_length_records
    if retained_records == 0:
        raise ValidationError(
            f"{path}: no records remain after applying N and length exclusions."
        )
    records = _records_from_abundances(abundance_by_sequence)
    return records, LibraryReadSummary(
        input_records=input_records,
        input_total_abundance=input_records,
        excluded_n_records=excluded_n_records,
        excluded_n_abundance=excluded_n_records,
        excluded_length_records=excluded_length_records,
        excluded_length_abundance=excluded_length_records,
        retained_records=retained_records,
        retained_unique_sequences=len(records),
        retained_abundance=retained_records,
    )


def read_collapsed_tag(
    path: Path,
    *,
    min_length: int | None = None,
    max_length: int | None = None,
) -> tuple[list[FastqRead], LibraryReadSummary]:
    """Read a two-column ``sequence<TAB>abundance`` table.

    A complete row is excluded when its sequence contains N. Other ambiguous
    symbols are errors, as are non-positive or non-integer abundance values.
    Duplicate canonical sequences are summed without encoding abundance in the
    Bowtie query identifier.
    """

    _validate_length_bounds(min_length, max_length)
    if path.suffix != ".tag":
        raise ValidationError(f"{path}: expected a collapsed .tag file.")
    if not path.is_file():
        raise ValidationError(f"Collapsed-tag library does not exist: {path}")

    abundance_by_sequence: dict[str, int] = {}
    input_records = 0
    input_total_abundance = 0
    excluded_n_records = 0
    excluded_n_abundance = 0
    excluded_length_records = 0
    excluded_length_abundance = 0
    retained_records = 0
    retained_abundance = 0
    with _open_text(path) as handle:
        for row_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line:
                raise ValidationError(f"{path}:{row_number}: blank tag rows are not allowed.")
            fields = line.split("\t")
            if len(fields) != 2:
                raise ValidationError(
                    f"{path}:{row_number}: expected exactly sequence<TAB>abundance."
                )
            raw_sequence, raw_abundance = fields
            try:
                abundance = int(raw_abundance.strip())
            except ValueError as exc:
                raise ValidationError(
                    f"{path}:{row_number}: abundance must be a positive integer."
                ) from exc
            if abundance <= 0:
                raise ValidationError(
                    f"{path}:{row_number}: abundance must be a positive integer."
                )
            input_records += 1
            input_total_abundance += abundance
            normalized = _normalized_input_sequence(raw_sequence)
            if "N" in normalized:
                excluded_n_records += 1
                excluded_n_abundance += abundance
                continue
            sequence = canonical_sequence(normalized, context=f"{path}:{row_number}")
            if _outside_length_bounds(sequence, min_length, max_length):
                excluded_length_records += 1
                excluded_length_abundance += abundance
                continue
            retained_records += 1
            retained_abundance += abundance
            abundance_by_sequence[sequence] = (
                abundance_by_sequence.get(sequence, 0) + abundance
            )
    if input_records == 0:
        raise ValidationError(f"{path}: no tag rows found.")
    if retained_records == 0:
        raise ValidationError(
            f"{path}: no rows remain after applying N and length exclusions."
        )
    records = _records_from_abundances(abundance_by_sequence)
    return records, LibraryReadSummary(
        input_records=input_records,
        input_total_abundance=input_total_abundance,
        excluded_n_records=excluded_n_records,
        excluded_n_abundance=excluded_n_abundance,
        excluded_length_records=excluded_length_records,
        excluded_length_abundance=excluded_length_abundance,
        retained_records=retained_records,
        retained_unique_sequences=len(records),
        retained_abundance=retained_abundance,
    )


def read_library_collapsed(
    path: Path,
    *,
    library_format: str,
    min_length: int | None = None,
    max_length: int | None = None,
) -> tuple[list[FastqRead], LibraryReadSummary]:
    """Read and collapse one explicitly declared library format."""

    if library_format == "fastq":
        return read_collapsed_fastq(
            path, min_length=min_length, max_length=max_length
        )
    if library_format == "tag":
        return read_collapsed_tag(path, min_length=min_length, max_length=max_length)
    raise ValidationError(
        f"Unsupported library format {library_format!r}; choose from: {', '.join(LIBRARY_SUFFIXES)}."
    )


def write_mapping_queries(records: Iterable[FastqRead], destination: Path) -> None:
    """Write canonical internal FASTA queries for Bowtie without abundance headers."""

    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("wt", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(f">{record.query_id}\n{record.sequence}\n")


def read_targets(path: Path) -> list[Sample]:
    """Read a tab-separated target sheet with explicit biological replication."""

    if not path.is_file():
        raise ValidationError(f"Target file does not exist: {path}")
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "condition", "biological_replicate"}
        fields = set(reader.fieldnames or [])
        missing = sorted(required - fields)
        if missing:
            raise ValidationError(
                f"{path}: missing required target columns: {', '.join(missing)}."
            )
        samples: list[Sample] = []
        seen: set[str] = set()
        for row_number, row in enumerate(reader, start=2):
            values = {key: (value or "").strip() for key, value in row.items()}
            sample_id = values["sample_id"]
            condition = values["condition"]
            replicate = values["biological_replicate"]
            if not sample_id or not condition or not replicate:
                raise ValidationError(
                    f"{path}:{row_number}: sample_id, condition, and biological_replicate "
                    "must all be non-empty."
                )
            if sample_id in seen:
                raise ValidationError(f"{path}:{row_number}: duplicated sample_id {sample_id!r}.")
            seen.add(sample_id)
            samples.append(
                Sample(
                    sample_id=sample_id,
                    condition=condition,
                    biological_replicate=replicate,
                    library_id=values.get("library_id") or sample_id,
                    batch=values.get("batch") or None,
                    technical_replicate=values.get("technical_replicate") or None,
                )
            )
    if not samples:
        raise ValidationError(f"{path}: no target rows found.")
    return samples


def find_processed_libraries(
    processed_libraries: Path, *, library_format: str = "fastq"
) -> dict[str, Path]:
    """Find individual libraries of one explicitly declared format."""

    if not processed_libraries.is_dir():
        raise ValidationError(f"Processed-libraries directory does not exist: {processed_libraries}")
    try:
        suffix = LIBRARY_SUFFIXES[library_format]
    except KeyError as exc:
        raise ValidationError(
            f"Unsupported library format {library_format!r}; choose from: {', '.join(LIBRARY_SUFFIXES)}."
        ) from exc
    found: dict[str, Path] = {}
    for candidate in sorted(processed_libraries.glob(f"*{suffix}")):
        found[library_stem(candidate, library_format=library_format)] = candidate
    if not found:
        raise ValidationError(
            f"{processed_libraries}: no {suffix} libraries found for --library-format {library_format}."
        )
    return found


def library_stem(path: Path, *, library_format: str = "fastq") -> str:
    """Return the basename without the suffix for the declared library format."""

    name = path.name
    try:
        suffix = LIBRARY_SUFFIXES[library_format]
    except KeyError as exc:
        raise ValidationError(
            f"Unsupported library format {library_format!r}; choose from: {', '.join(LIBRARY_SUFFIXES)}."
        ) from exc
    if name.endswith(suffix):
        return name[: -len(suffix)]
    raise ValidationError(f"{path}: expected a {suffix} suffix.")


def resolve_target_libraries(
    samples: Iterable[Sample],
    processed_libraries: Path,
    *,
    library_format: str = "fastq",
) -> dict[str, Path]:
    """Resolve every target row to exactly one library of the declared format."""

    available = find_processed_libraries(
        processed_libraries, library_format=library_format
    )
    suffix = LIBRARY_SUFFIXES[library_format]
    resolved: dict[str, Path] = {}
    for sample in samples:
        path = available.get(sample.library_id)
        if path is None:
            suggestions = ", ".join(sorted(available)[:12])
            raise ValidationError(
                f"Target sample {sample.sample_id!r} requests library_id {sample.library_id!r}, "
                f"but no matching <library>{suffix} exists in {processed_libraries}. "
                f"Available examples: {suggestions}"
            )
        resolved[sample.sample_id] = path
    return resolved


def read_reference_fasta(path: Path) -> list[tuple[str, str]]:
    """Read a FASTA reference and return ``(feature_id, canonical_sequence)`` records."""

    if not path.is_file():
        raise ValidationError(f"Reference FASTA does not exist: {path}")
    records: list[tuple[str, str]] = []
    seen_ids: set[str] = set()
    header: str | None = None
    sequence_lines: list[str] = []

    def emit(current_header: str, lines: list[str]) -> None:
        feature_id = current_header[1:].strip().split(maxsplit=1)[0] if current_header[1:].strip() else ""
        if not feature_id:
            raise ValidationError(f"{path}: FASTA header {current_header!r} has no feature ID.")
        if feature_id in seen_ids:
            raise ValidationError(f"{path}: duplicate FASTA feature ID {feature_id!r}.")
        if not lines:
            raise ValidationError(f"{path}: feature {feature_id!r} has no sequence.")
        seen_ids.add(feature_id)
        records.append(
            (feature_id, canonical_sequence("".join(lines), context=f"{path} {feature_id}"))
        )

    with _open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line:
                raise ValidationError(f"{path}:{line_number}: blank lines are not allowed in FASTA.")
            if line.startswith(">"):
                if header is not None:
                    emit(header, sequence_lines)
                header = line
                sequence_lines = []
            else:
                if header is None:
                    raise ValidationError(f"{path}:{line_number}: sequence before FASTA header.")
                sequence_lines.append(line)
    if header is not None:
        emit(header, sequence_lines)
    if not records:
        raise ValidationError(f"{path}: no FASTA records found.")
    return records


def parse_reference_spec(value: str) -> tuple[str, Path]:
    """Parse ``CLASS=PATH`` command-line reference declarations."""

    if "=" not in value:
        raise ValidationError(
            f"Invalid --reference value {value!r}; expected CLASS=/path/to/reference.fasta."
        )
    srna_class, raw_path = value.split("=", maxsplit=1)
    srna_class = srna_class.strip()
    if not SAFE_CLASS.fullmatch(srna_class):
        raise ValidationError(
            f"Invalid reference class {srna_class!r}; use letters, digits, '_' or '-'."
        )
    if not raw_path.strip():
        raise ValidationError(f"Invalid --reference value {value!r}; path is empty.")
    return srna_class, Path(raw_path).expanduser().resolve()


def sha256_file(path: Path) -> str:
    """Return a streaming SHA-256 checksum without loading an input into memory."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, object]]) -> None:
    """Write a stable tab-separated table with a header, creating parents as needed."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})
