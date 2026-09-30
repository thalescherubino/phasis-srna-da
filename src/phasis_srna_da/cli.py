"""Command-line entry point for conservative small-RNA counting."""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from . import __version__
from .catalog import (
    Catalog,
    build_catalog,
    catalog_decision_rows,
    catalog_record_rows,
    write_atomic_reference,
)
from .counts import CountAccumulator, Hierarchy, read_hierarchy, write_gzip_tsv
from .errors import ExternalToolError, ValidationError
from .io import (
    LibraryReadSummary,
    parse_reference_spec,
    read_library_collapsed,
    read_targets,
    resolve_target_libraries,
    sha256_file,
    write_mapping_queries,
    write_tsv,
)
from .mapping import (
    assignments_from_sam,
    build_bowtie_index,
    executable_version,
    run_bowtie_exact,
)
from .models import Sample
from .provenance import (
    base_run_metadata,
    prepare_output_directory,
    utc_timestamp,
    write_input_manifest,
    write_json,
)


LEDGER_COLUMNS = [
    "sample_id",
    "query_id",
    "sequence",
    "abundance",
    "atomic_ids",
    "atomic_features",
    "primary_assignment_type",
    "primary_unit_id",
]

SUPPORTED_REFERENCE_ORIENTATIONS = ("molecule",)
SUPPORTED_NORMALIZATION_POLICIES = ("total-qc-depth",)
SUPPORTED_PREFILTER_POLICIES = (
    "all-replicates-in-one-condition",
    "any-samples",
)
SUPPORTED_LIBRARY_FORMATS = ("fastq", "tag")
AMBIGUOUS_N_POLICY = "drop_complete_library_record"


@dataclass(frozen=True)
class PreparedInputs:
    samples: tuple[Sample, ...]
    library_paths: dict[str, Path]
    catalog: Catalog
    hierarchy: Hierarchy | None
    assembly_id: str
    annotation_id: str
    library_format: str
    min_length: int
    max_length: int
    reference_orientation: str
    library_summaries: dict[str, LibraryReadSummary]
    warnings: tuple[str, ...]


def _path(value: str) -> Path:
    return Path(value).expanduser().resolve()


def _identity(value: str) -> str:
    identity = value.strip()
    if not identity:
        raise argparse.ArgumentTypeError("identity must not be blank")
    return identity


def _positive_int(value: str) -> int:
    try:
        number = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("value must be a positive integer") from exc
    if number < 1:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return number


def _contrast(value: str) -> tuple[str, str]:
    """Parse one explicit ``NUMERATOR,DENOMINATOR`` contrast declaration."""

    fields = [field.strip() for field in value.split(",")]
    if len(fields) != 2 or not all(fields):
        raise argparse.ArgumentTypeError(
            "contrast must be NUMERATOR_CONDITION,DENOMINATOR_CONDITION"
        )
    numerator, denominator = fields
    if numerator == denominator:
        raise argparse.ArgumentTypeError("contrast conditions must differ")
    return numerator, denominator


def _add_shared_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--processed-libraries",
        type=_path,
        required=True,
        metavar="DIR",
        help=(
            "Directory containing one <library_id>.fastq.gz or <library_id>.tag file "
            "per library, as selected by --library-format."
        ),
    )
    parser.add_argument(
        "--library-format",
        required=True,
        choices=SUPPORTED_LIBRARY_FORMATS,
        help=(
            "Explicit input format: 'fastq' for .fastq.gz records or 'tag' for "
            "two-column sequence<TAB>abundance files. Records containing N are dropped."
        ),
    )
    parser.add_argument(
        "--min-length",
        required=True,
        type=_positive_int,
        metavar="NT",
        help="Minimum retained library sequence length, declared explicitly in nucleotides.",
    )
    parser.add_argument(
        "--max-length",
        required=True,
        type=_positive_int,
        metavar="NT",
        help="Maximum retained library sequence length, declared explicitly in nucleotides.",
    )
    parser.add_argument(
        "--targets",
        type=_path,
        required=True,
        metavar="TSV",
        help=(
            "Tab-separated targets with sample_id, condition, biological_replicate; "
            "optional library_id, batch, technical_replicate."
        ),
    )
    parser.add_argument(
        "--assembly-id",
        required=True,
        type=_identity,
        metavar="ID",
        help="Genome assembly identity used to generate the supplied Phasis products.",
    )
    parser.add_argument(
        "--annotation-id",
        required=True,
        type=_identity,
        metavar="ID",
        help="Stable identity of the frozen sRNA annotation/reference collection.",
    )
    parser.add_argument(
        "--reference",
        action="append",
        required=True,
        metavar="CLASS=FASTA",
        help="Repeat for each supplied class (for example miRNA=/refs/mature.fa).",
    )
    parser.add_argument(
        "--phasiRNAs",
        "--phasi-catalog",
        dest="phasi_catalog",
        action="append",
        required=True,
        type=_path,
        metavar="TSV",
        help=(
            "Required Phasis *_phasiRNAs.tsv product file. Repeat for each explicitly "
            "selected fixed 21- or 24-nt product file."
        ),
    )
    parser.add_argument(
        "--reference-orientation",
        required=True,
        choices=SUPPORTED_REFERENCE_ORIENTATIONS,
        help=(
            "Required input contract. 'molecule' means supplied FASTA sequences are "
            "written 5'-to-3' in the same molecular orientation as the reads; mapping "
            "is direct only and does not search reverse complements."
        ),
    )
    parser.add_argument(
        "--hierarchy",
        type=_path,
        metavar="TSV",
        help="Optional class/feature_id/parent_level/parent_id hierarchy metadata.",
    )
    parser.add_argument(
        "--allow-catalog-union",
        action="store_true",
        help=(
            "Allow multiple explicitly supplied Phasis catalogs for one phase. "
            "Use only for a documented, curated fixed union."
        ),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="phasis-srna-da",
        description=(
            "Re-count mandatory fixed Phasis product files together with explicit "
            "small-RNA reference FASTAs in explicitly declared FASTQ or collapsed-tag libraries. "
            "No de novo locus discovery is performed."
        ),
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate = subparsers.add_parser("validate", help="Validate all inputs without Bowtie mapping.")
    _add_shared_arguments(validate)
    validate.add_argument(
        "--json",
        action="store_true",
        help="Emit the validation summary as JSON rather than human-readable text.",
    )

    run = subparsers.add_parser("run", help="Build an index, re-count libraries, and write DA-ready matrices.")
    _add_shared_arguments(run)
    run.add_argument("--outdir", type=_path, required=True, metavar="DIR")
    run.add_argument("--threads", type=int, default=1, help="Bowtie 1 threads (default: 1).")
    run.add_argument("--bowtie", default="bowtie", help="Bowtie 1 executable or PATH name.")
    run.add_argument(
        "--bowtie-build", default="bowtie-build", help="Bowtie-build executable or PATH name."
    )
    run.add_argument(
        "--keep-alignments",
        action="store_true",
        help="Keep raw SAM files under outdir/work/alignments after audited parsing.",
    )
    run.add_argument(
        "--run-da",
        action="store_true",
        help=(
            "Run optional PyDESeq2 models using one fixed retained-library depth factor "
            "per sample across every count projection."
        ),
    )
    run.add_argument(
        "--normalization",
        choices=SUPPORTED_NORMALIZATION_POLICIES,
        help=(
            "Required with --run-da. 'total-qc-depth' uses one size factor derived from "
            "the retained post-N-exclusion abundance of each biological sample and "
            "shares it across all count projections."
        ),
    )
    run.add_argument(
        "--contrast",
        action="append",
        type=_contrast,
        metavar="NUMERATOR,DENOMINATOR",
        help=(
            "Required with --run-da; repeat for every intended comparison. Conditions "
            "must exactly match --targets. No undeclared pairwise contrasts are run."
        ),
    )
    run.add_argument("--min-count", type=int, default=10, help="DA prefilter raw count (default: 10).")
    run.add_argument(
        "--prefilter-policy", choices=SUPPORTED_PREFILTER_POLICIES,
        default="all-replicates-in-one-condition",
        help=(
            "DA low-count rule (default: all-replicates-in-one-condition, which requires "
            "--min-count in every biological replicate of at least one condition)."
        ),
    )
    run.add_argument("--fdr", type=float, default=0.05, help="DA FDR alpha (default: 0.05).")

    da = subparsers.add_parser(
        "da", help="Run explicit contrasts from a completed Phasis-backed count run without remapping."
    )
    da.add_argument(
        "--count-run", type=_path, required=True, metavar="DIR",
        help="Completed 'run' directory; inherits its targets, mandatory Phasis provenance, and retained depths.",
    )
    da.add_argument("--outdir", type=_path, required=True, metavar="DIR")
    da.add_argument(
        "--normalization", choices=SUPPORTED_NORMALIZATION_POLICIES, required=True,
        help="Use total retained library abundance within the source run's declared length range.",
    )
    da.add_argument(
        "--contrast", action="append", type=_contrast, required=True,
        metavar="NUMERATOR,DENOMINATOR",
        help="Repeat for each comparison; condition names must match the count run's targets.",
    )
    da.add_argument("--threads", type=_positive_int, default=1)
    da.add_argument("--min-count", type=_positive_int, default=10)
    da.add_argument(
        "--prefilter-policy", choices=SUPPORTED_PREFILTER_POLICIES,
        default="all-replicates-in-one-condition",
    )
    da.add_argument("--fdr", type=float, default=0.05)
    return parser


def _validate_run_policy(args: argparse.Namespace) -> None:
    """Require an explicit DA normalization choice without implying DA for count-only runs."""

    if args.run_da and args.normalization is None:
        raise ValidationError(
            "--normalization is required with --run-da; currently supported: total-qc-depth."
        )
    if not args.run_da and args.normalization is not None:
        raise ValidationError("--normalization is only valid with --run-da.")
    contrasts = getattr(args, "contrast", None) or []
    if args.run_da and not contrasts:
        raise ValidationError(
            "At least one --contrast is required with --run-da; use "
            "NUMERATOR_CONDITION,DENOMINATOR_CONDITION."
        )
    if not args.run_da and contrasts:
        raise ValidationError("--contrast is only valid with --run-da.")


def _validate_declared_contrasts(
    contrasts: list[tuple[str, str]], conditions: set[str]
) -> None:
    """Require unique, directed contrasts between conditions in the target sheet."""

    seen: set[tuple[str, str]] = set()
    for numerator, denominator in contrasts:
        missing = sorted({numerator, denominator} - conditions)
        if missing:
            raise ValidationError(
                "Contrast condition(s) not present in --targets: " + ", ".join(missing) + "."
            )
        directed = (numerator, denominator)
        if directed in seen:
            raise ValidationError(
                f"Duplicate contrast declaration: {numerator},{denominator}."
            )
        if (denominator, numerator) in seen:
            raise ValidationError(
                "Both directions of one contrast were declared; choose one numerator: "
                f"{numerator},{denominator}."
            )
        seen.add(directed)


def _target_warnings(samples: tuple[Sample, ...]) -> list[str]:
    warnings: list[str] = []
    condition_counts = Counter(sample.condition for sample in samples)
    if len(condition_counts) < 2:
        warnings.append("Only one condition is present: matrices can be built but DA contrasts are unavailable.")
    for condition, count in sorted(condition_counts.items()):
        if count < 2:
            warnings.append(
                f"Condition {condition!r} has {count} library; inferential DA needs at least 2 independent replicates."
            )
        elif count < 3:
            warnings.append(
                f"Condition {condition!r} has {count} libraries; 3 biological replicates are recommended for publication."
            )
    seen_biological: set[tuple[str, str]] = set()
    for sample in samples:
        key = (sample.condition, sample.biological_replicate)
        if key in seen_biological:
            warnings.append(
                "Repeated condition/biological_replicate values detected. Confirm these are technical "
                "replicates and combine them before inferential DA."
            )
            break
        seen_biological.add(key)
    if any(sample.batch for sample in samples) and not all(sample.batch for sample in samples):
        warnings.append("Batch is partially specified; no batch-adjusted model will be used until it is complete.")
    return warnings


def prepare_inputs(args: argparse.Namespace, *, validate_libraries: bool) -> PreparedInputs:
    """Validate fixed input contracts before writing any mapping result."""

    if getattr(args, "threads", 1) < 1:
        raise ValidationError("--threads must be at least 1.")
    if args.min_length > args.max_length:
        raise ValidationError("--min-length must not exceed --max-length.")
    samples = tuple(read_targets(args.targets))
    library_paths = resolve_target_libraries(
        samples,
        args.processed_libraries,
        library_format=args.library_format,
    )
    catalog = build_catalog(
        args.reference,
        list(args.phasi_catalog),
        allow_catalog_union=args.allow_catalog_union,
    )
    hierarchy = read_hierarchy(args.hierarchy) if args.hierarchy else None
    library_summaries: dict[str, LibraryReadSummary] = {}
    if validate_libraries:
        for sample in samples:
            _, library_summaries[sample.sample_id] = read_library_collapsed(
                library_paths[sample.sample_id],
                library_format=args.library_format,
                min_length=args.min_length,
                max_length=args.max_length,
            )
    return PreparedInputs(
        samples=samples,
        library_paths=library_paths,
        catalog=catalog,
        hierarchy=hierarchy,
        assembly_id=args.assembly_id,
        annotation_id=args.annotation_id,
        library_format=args.library_format,
        min_length=args.min_length,
        max_length=args.max_length,
        reference_orientation=args.reference_orientation,
        library_summaries=library_summaries,
        warnings=tuple(_target_warnings(samples)),
    )


def _validation_summary(prepared: PreparedInputs) -> dict[str, object]:
    selected = sum(decision.selected for decision in prepared.catalog.phasi_decisions)
    classes = Counter(record.srna_class for record in prepared.catalog.records)
    summary: dict[str, object] = {
        "valid": True,
        "samples": len(prepared.samples),
        "conditions": sorted({sample.condition for sample in prepared.samples}),
        "assembly_id": prepared.assembly_id,
        "annotation_id": prepared.annotation_id,
        "library_format": prepared.library_format,
        "library_min_length": prepared.min_length,
        "library_max_length": prepared.max_length,
        "ambiguous_N_policy": AMBIGUOUS_N_POLICY,
        "atomic_reference_records": len(prepared.catalog.records),
        "atomic_reference_records_by_class": dict(sorted(classes.items())),
        "phasi_product_files": len(
            {
                decision.association.source_path
                for decision in prepared.catalog.phasi_decisions
            }
        ),
        "phasi_catalog_rows_selected": selected,
        "phasi_catalog_rows_seen": len(prepared.catalog.phasi_decisions),
        "reference_orientation": prepared.reference_orientation,
        "mapping_orientation": "direct_only_no_reverse_complement",
        "warnings": list(prepared.warnings),
    }
    if prepared.library_summaries:
        summaries = tuple(prepared.library_summaries.values())
        summary["library_input_records"] = sum(item.input_records for item in summaries)
        summary["library_input_total_abundance"] = sum(
            item.input_total_abundance for item in summaries
        )
        summary["library_excluded_N_records"] = sum(
            item.excluded_n_records for item in summaries
        )
        summary["library_excluded_N_abundance"] = sum(
            item.excluded_n_abundance for item in summaries
        )
        summary["library_excluded_length_records"] = sum(
            item.excluded_length_records for item in summaries
        )
        summary["library_excluded_length_abundance"] = sum(
            item.excluded_length_abundance for item in summaries
        )
        summary["library_retained_records"] = sum(
            item.retained_records for item in summaries
        )
        summary["library_retained_unique_sequences"] = sum(
            item.retained_unique_sequences for item in summaries
        )
        summary["library_retained_abundance"] = sum(
            item.retained_abundance for item in summaries
        )
    return summary


def command_validate(args: argparse.Namespace) -> int:
    prepared = prepare_inputs(args, validate_libraries=True)
    summary = _validation_summary(prepared)
    if args.json:
        print(json.dumps(summary, indent=2, sort_keys=True))
    else:
        print("Inputs are valid.")
        print(f"  samples: {summary['samples']}")
        print(f"  conditions: {', '.join(summary['conditions'])}")
        print(f"  assembly: {summary['assembly_id']}")
        print(f"  annotation: {summary['annotation_id']}")
        print(f"  library format: {summary['library_format']}")
        print(
            "  retained length range: "
            f"{summary['library_min_length']}--{summary['library_max_length']} nt"
        )
        print(f"  N policy: {summary['ambiguous_N_policy']}")
        print(f"  atomic reference records: {summary['atomic_reference_records']}")
        print(f"  required Phasis product files: {summary['phasi_product_files']}")
        print(f"  selected Phasis catalog rows: {summary['phasi_catalog_rows_selected']}")
        print(f"  reference orientation: {summary['reference_orientation']}")
        for warning in summary["warnings"]:
            print(f"  warning: {warning}")
    return 0


def _input_manifest_entries(args: argparse.Namespace, prepared: PreparedInputs) -> list[tuple[str, Path, str]]:
    entries: list[tuple[str, Path, str]] = [("targets", args.targets, "target/design sheet")]
    for sample in prepared.samples:
        role = (
            "quality_controlled_fastq"
            if prepared.library_format == "fastq"
            else "collapsed_tag_table"
        )
        entries.append(
            (
                role,
                prepared.library_paths[sample.sample_id],
                (
                    f"sample_id={sample.sample_id}; format={prepared.library_format}; "
                    "drop complete records containing N; "
                    f"retain {prepared.min_length}--{prepared.max_length} nt; "
                    "preserve retained abundance"
                ),
            )
        )
    for spec in args.reference:
        srna_class, path = parse_reference_spec(spec)
        entries.append(
            (
                "reference_fasta",
                path,
                f"class={srna_class}; orientation={prepared.reference_orientation}",
            )
        )
    for path in args.phasi_catalog:
        entries.append(
            (
                "phasis_product_catalog",
                path,
                "fixed user-supplied Phasis product catalog; tag_seq is molecule-oriented; no discovery",
            )
        )
    if args.hierarchy:
        entries.append(("optional_hierarchy", args.hierarchy, "declared parent relationships"))
    return entries


def _write_catalog_outputs(outdir: Path, catalog: Catalog) -> None:
    write_tsv(
        outdir / "catalogs" / "atomic_reference_features.tsv",
        [
            "atomic_id",
            "srna_class",
            "feature_id",
            "sequence",
            "source_path",
            "source_kind",
            "phase",
            "locus_id",
            "tag_sequence",
            "strand",
            "observed_pos",
        ],
        catalog_record_rows(catalog.records),
    )
    decision_fields = [
        "selected",
        "decision",
        "phase",
        "identifier",
        "tag_seq",
        "register_class",
        "window_unit_id",
        "window_unit_role",
        "window_unit_rank",
        "window_unit_shift_nt",
        "strand",
        "expected_register_pos",
        "observed_pos",
        "cID",
        "alib",
        "abun",
        "hits",
        "source_path",
        "source_row",
        "source_row_sha256",
        "physical_product_key",
        "atomic_id",
    ]
    write_tsv(
        outdir / "catalogs" / "phasi_product_catalog_decisions.tsv",
        decision_fields,
        catalog_decision_rows(catalog),
    )


def _write_sample_metadata(
    outdir: Path,
    samples: tuple[Sample, ...],
    library_paths: dict[str, Path],
    library_format: str,
) -> None:
    write_tsv(
        outdir / "provenance" / "resolved_targets.tsv",
        [
            "sample_id",
            "condition",
            "biological_replicate",
            "library_id",
            "batch",
            "technical_replicate",
            "library_format",
            "library_path",
        ],
        [
            {
                "sample_id": sample.sample_id,
                "condition": sample.condition,
                "biological_replicate": sample.biological_replicate,
                "library_id": sample.library_id,
                "batch": sample.batch or "",
                "technical_replicate": sample.technical_replicate or "",
                "library_format": library_format,
                "library_path": str(library_paths[sample.sample_id]),
            }
            for sample in samples
        ],
    )


def _write_count_outputs(outdir: Path, accumulator: CountAccumulator) -> None:
    for projection_name in sorted(accumulator.projections):
        fields, rows = accumulator.count_matrix_rows(projection_name)
        write_tsv(outdir / "counts" / f"{projection_name}.raw.tsv", fields, rows)
    fields, rows = accumulator.equivalence_rows()
    write_tsv(outdir / "mapping" / "equivalence_group_members.tsv", fields, rows)
    fields, rows = accumulator.conservation_rows()
    write_tsv(outdir / "mapping" / "conservation_by_sample.tsv", fields, rows)


def command_run(args: argparse.Namespace) -> int:
    if args.min_count < 1:
        raise ValidationError("--min-count must be at least 1.")
    if not 0 < args.fdr < 1:
        raise ValidationError("--fdr must be between 0 and 1.")
    _validate_run_policy(args)
    prepared = prepare_inputs(args, validate_libraries=False)
    if args.run_da:
        _validate_declared_contrasts(
            args.contrast,
            {sample.condition for sample in prepared.samples},
        )
    prepare_output_directory(args.outdir)
    metadata = base_run_metadata(
        command_line=getattr(args, "_command_line", list(sys.argv))
    )
    metadata["validation"] = _validation_summary(prepared)
    metadata["allow_catalog_union"] = args.allow_catalog_union
    metadata["assembly_id"] = prepared.assembly_id
    metadata["annotation_id"] = prepared.annotation_id
    metadata["library_format"] = prepared.library_format
    metadata["library_min_length"] = prepared.min_length
    metadata["library_max_length"] = prepared.max_length
    metadata["ambiguous_N_policy"] = AMBIGUOUS_N_POLICY
    metadata["reference_orientation"] = prepared.reference_orientation
    metadata["mapping_orientation"] = "direct_only_no_reverse_complement"
    metadata["da_normalization_policy"] = args.normalization if args.run_da else None
    metadata["da_prefilter_policy"] = args.prefilter_policy if args.run_da else None
    metadata["da_contrasts"] = (
        [
            {"numerator": numerator, "denominator": denominator}
            for numerator, denominator in args.contrast
        ]
        if args.run_da
        else []
    )
    metadata["status"] = "running"
    write_json(args.outdir / "provenance" / "run_metadata.json", metadata)
    write_input_manifest(args.outdir, _input_manifest_entries(args, prepared))
    _write_sample_metadata(
        args.outdir,
        prepared.samples,
        prepared.library_paths,
        prepared.library_format,
    )
    _write_catalog_outputs(args.outdir, prepared.catalog)
    write_tsv(
        args.outdir / "provenance" / "warnings.tsv",
        ["warning"],
        [{"warning": warning} for warning in prepared.warnings],
    )

    reference_fasta = args.outdir / "catalogs" / "atomic_reference.fa"
    write_atomic_reference(prepared.catalog.records, reference_fasta)
    bowtie_build_log = args.outdir / "provenance" / "bowtie-build.log"
    index_prefix = args.outdir / "work" / "bowtie_index" / "atomic_reference"
    index_run = build_bowtie_index(
        reference_fasta,
        index_prefix,
        bowtie_build=args.bowtie_build,
        stderr_path=bowtie_build_log,
    )
    accumulator = CountAccumulator(
        (sample.sample_id for sample in prepared.samples),
        prepared.catalog.records,
        hierarchy=prepared.hierarchy,
    )
    library_rows: list[dict[str, object]] = []
    mapping_commands: list[dict[str, str]] = [
        {
            "step": "bowtie-build",
            "command": index_run.rendered_command,
            "stderr_path": index_run.stderr_path,
        }
    ]
    for sample in prepared.samples:
        source_path = prepared.library_paths[sample.sample_id]
        print(f"[{sample.sample_id}] Reading and filtering library", file=sys.stderr, flush=True)
        records, library_summary = read_library_collapsed(
            source_path,
            library_format=prepared.library_format,
            min_length=prepared.min_length,
            max_length=prepared.max_length,
        )
        query_path = args.outdir / "work" / "queries" / f"{sample.sample_id}.fasta"
        write_mapping_queries(records, query_path)
        sam_path = args.outdir / "work" / "alignments" / f"{sample.sample_id}.sam"
        bowtie_log = args.outdir / "provenance" / "bowtie_logs" / f"{sample.sample_id}.log"
        print(
            f"[{sample.sample_id}] Mapping {len(records):,} sequences representing "
            f"{library_summary.retained_abundance:,} retained reads",
            file=sys.stderr,
            flush=True,
        )
        mapping_run = run_bowtie_exact(
            index_prefix,
            query_path,
            sam_path,
            threads=args.threads,
            bowtie=args.bowtie,
            stderr_path=bowtie_log,
        )
        mapping_commands.append(
            {
                "step": f"bowtie:{sample.sample_id}",
                "command": mapping_run.rendered_command,
                "stderr_path": mapping_run.stderr_path,
            }
        )
        assignments = assignments_from_sam(
            sample.sample_id, records, sam_path, prepared.catalog.records
        )
        write_gzip_tsv(
            args.outdir / "mapping" / "assignments" / f"{sample.sample_id}.tsv.gz",
            LEDGER_COLUMNS,
            accumulator.add_assignments(assignments),
        )
        primary_summary = accumulator.projections["global_atomic"].conservation[sample.sample_id]
        if primary_summary["input_abundance"] != library_summary.retained_abundance:
            raise ValidationError(f"{sample.sample_id}: mapped input abundance differs from retained library abundance.")
        library_rows.append(
            {
                "sample_id": sample.sample_id,
                "condition": sample.condition,
                "biological_replicate": sample.biological_replicate,
                "library_path": str(source_path),
                "library_format": prepared.library_format,
                "input_records": library_summary.input_records,
                "input_total_abundance": library_summary.input_total_abundance,
                "excluded_N_records": library_summary.excluded_n_records,
                "excluded_N_abundance": library_summary.excluded_n_abundance,
                "excluded_length_records": library_summary.excluded_length_records,
                "excluded_length_abundance": library_summary.excluded_length_abundance,
                "retained_records": library_summary.retained_records,
                "retained_unique_sequences": library_summary.retained_unique_sequences,
                "retained_abundance": library_summary.retained_abundance,
                "normalization_depth": library_summary.retained_abundance,
                "primary_unique_assigned_abundance": primary_summary["unique_assigned_abundance"],
                "primary_equivalence_assigned_abundance": primary_summary[
                    "equivalence_assigned_abundance"
                ],
                "primary_unassigned_abundance": primary_summary["unassigned_abundance"],
            }
        )
        if not args.keep_alignments:
            sam_path.unlink()
        write_tsv(
            args.outdir / "qc" / "library_summary.tsv",
            list(library_rows[0]),
            library_rows,
        )
        print(
            f"[{sample.sample_id}] Counted; retained abundance conserved",
            file=sys.stderr,
            flush=True,
        )
        del records, assignments

    _write_count_outputs(args.outdir, accumulator)
    write_tsv(
        args.outdir / "qc" / "library_summary.tsv",
        [
            "sample_id",
            "condition",
            "biological_replicate",
            "library_path",
            "library_format",
            "input_records",
            "input_total_abundance",
            "excluded_N_records",
            "excluded_N_abundance",
            "excluded_length_records",
            "excluded_length_abundance",
            "retained_records",
            "retained_unique_sequences",
            "retained_abundance",
            "normalization_depth",
            "primary_unique_assigned_abundance",
            "primary_equivalence_assigned_abundance",
            "primary_unassigned_abundance",
        ],
        library_rows,
    )
    write_tsv(
        args.outdir / "provenance" / "commands.tsv",
        ["step", "command", "stderr_path"],
        mapping_commands,
    )
    write_json(
        args.outdir / "provenance" / "software_versions.json",
        {
            "bowtie": executable_version(args.bowtie),
            "bowtie_build": executable_version(args.bowtie_build),
            "phasis_srna_da": metadata["tool_version"],
        },
    )
    if args.run_da:
        from .da import run_differential_abundance

        run_differential_abundance(
            count_directory=args.outdir / "counts",
            samples=prepared.samples,
            outdir=args.outdir / "da",
            min_count=args.min_count,
            alpha=args.fdr,
            threads=args.threads,
            library_summary_path=args.outdir / "qc" / "library_summary.tsv",
            contrasts=tuple(args.contrast),
            prefilter_policy=args.prefilter_policy,
        )
    metadata["library_qc_totals"] = {
        "input_records": sum(int(row["input_records"]) for row in library_rows),
        "input_total_abundance": sum(
            int(row["input_total_abundance"]) for row in library_rows
        ),
        "excluded_N_records": sum(
            int(row["excluded_N_records"]) for row in library_rows
        ),
        "excluded_N_abundance": sum(
            int(row["excluded_N_abundance"]) for row in library_rows
        ),
        "excluded_length_records": sum(
            int(row["excluded_length_records"]) for row in library_rows
        ),
        "excluded_length_abundance": sum(
            int(row["excluded_length_abundance"]) for row in library_rows
        ),
        "retained_abundance": sum(
            int(row["retained_abundance"]) for row in library_rows
        ),
    }
    metadata["status"] = "complete"
    metadata["completed_utc"] = utc_timestamp()
    metadata["atomic_reference_sha256"] = sha256_file(reference_fasta)
    write_json(args.outdir / "provenance" / "run_metadata.json", metadata)
    print(f"Completed count-conserving re-counting: {args.outdir}")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Run the CLI and convert expected input/tool failures into concise messages."""

    parser = build_parser()
    args = parser.parse_args(argv)
    args._command_line = (
        [parser.prog, *argv] if argv is not None else list(sys.argv)
    )
    try:
        if args.command == "validate":
            return command_validate(args)
        if args.command == "run":
            return command_run(args)
        if args.command == "da":
            from .resume_da import command_da

            return command_da(args)
        raise AssertionError(f"Unexpected command {args.command!r}")
    except (ValidationError, ExternalToolError, FileExistsError) as exc:
        parser.error(str(exc))
    return 2
