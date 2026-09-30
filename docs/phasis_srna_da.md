# `phasis-srna-da` input and analysis contract

## Purpose

`phasis-srna-da` is downstream of **Phasis**. It does not alter Phasis, rerun discovery, or infer missing sRNA-producing loci. Its count input is quality-controlled FASTQ or a two-column collapsed-tag table, selected explicitly. Phasis product catalogs are mandatory reference inputs; their recorded discovery abundance never substitutes for counts from the analysis libraries.

It creates integer count matrices and audit files suitable for Python differential-abundance models. The command-line tool is generic. Start with [software setup on Hive](phasis_srna_da_setup.md), then follow the [analysis tutorial](ms28_cli_tutorial.md) with your own dataset or the separate B73 WT/ms28 worked example.

## Installation

Use a separate analysis environment, never the environment used to run Phasis. On Hive, install packages and run analyses in a Slurm compute allocation, with environments and caches on group storage. Reuse an existing working environment when available:

```bash
git clone --branch v0.2.0.dev0 --single-branch https://github.com/thalescherubino/phasis-srna-da.git
cd phasis-srna-da
python -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install '.[da]'
.venv/bin/python -m phasis_srna_da --version
```

Bowtie **1** and `bowtie-build` must be available on `PATH`, or passed explicitly with `--bowtie` and `--bowtie-build`.

## Required inputs

### Target sheet

`targets.tsv` is tab-separated and requires:

```text
sample_id    condition    biological_replicate
Epi_1    Epi    1
Epi_2    Epi    2
```

Optional columns are `library_id`, `batch`, and `technical_replicate`. `library_id` is useful when a biological sample name differs from the FASTQ basename. A duplicate biological replicate is only permitted for explicitly labelled technical libraries, which are combined before DA.

At least two independent biological replicates per condition are required for inferential DA; three are the practical minimum for publication.

### Quality-controlled FASTQ libraries

With `--library-format fastq`, pass one directory containing one quality-controlled `.fastq.gz` file per target:

```text
processed_libraries/
  Epi_1.fastq.gz
  Epi_2.fastq.gz
  ...
```

The basename must match the `library_id` in `targets.tsv`; for example, `library_id=Epi_1` resolves to `Epi_1.fastq.gz`. Each library must use standard four-line FASTQ records:

```text
@read_identifier
AGAAGAGAGAGAGTACAGCCT
+
IIIIIIIIIIIIIIIIIIIII
```

Every retained FASTQ record counts once. The program collapses identical sequences internally while preserving their total abundance and creates `read_N` query IDs without abundance values. It rejects malformed FASTQ and sequence/quality length mismatches. Complete library records containing N are dropped; other ambiguous bases are rejected. It does not remove duplicate-read abundance, infer counts from headers, or perform adapter/quality trimming or contaminant removal.

With `--library-format tag`, each `<library_id>.tag` file must contain exactly `sequence<TAB>positive_integer_abundance`, without a header. Duplicate canonical sequences are summed. An N-containing row loses its entire abundance. Both formats require explicit inclusive `--min-length` and `--max-length` bounds; N exclusion precedes length exclusion. Case is normalized and U becomes T. References themselves must contain unambiguous sequences.

Preserve complete adapter-clean libraries and keep discovery inputs separate from DA inputs. For example, an 18--50-nt DA workflow test can include structural-RNA fragments while exact Phasis products remain 21 or 24 nt. Reference matching alone does not validate a longer fragment's biological role. Changing the length range also changes the retained-depth normalization denominator.

Use the same quality-control workflow for every sample, and use files in which structural RNA reads were retained. A `nocontam` library cannot support differential accumulation of rRNA or another class removed before counting.

### Class reference FASTAs

Supply each class with repeated `--reference CLASS=FASTA` options. The first whitespace-delimited token in a FASTA header is the feature ID; IDs must be unique within a class. References can contain full tRNAs, so a short exact sequence is allowed to align within a longer reference record.

The scope includes user-supplied miRNAs, tRNAs, snRNAs, snoRNAs, rRNAs, separately labelled candidate classes, and individual Phasis products. A generic FASTA match can be a substring of a supplied feature and should be described as reference-associated signal. Exact mature-product identity is enforced specifically for Phasis products. No de novo siRNA/repeat/TE discovery occurs.

Declare `--assembly-id`, `--annotation-id`, and `--reference-orientation molecule`. The orientation flag asserts that the supplied sequences are already written in the same 5'-to-3' molecular orientation as the reads; it does not convert genomic-strand exports. Use an audited, fixed molecule-oriented Phasis catalog where necessary. Do not change a frozen reference during a DA test.

### Fixed Phasis product catalogs

Supply one or more explicitly selected `*_phasiRNAs.tsv` files with `--phasi-catalog`. They define a fixed 21- and/or 24-nt product universe; their `abun`, `alib`, `cID`, and `hits` values are provenance only, never DA counts.

`--phasiRNAs` is the preferred spelling; `--phasi-catalog` is an alias. Multiple files for the same phase require explicit `--allow-catalog-union` and a documented fixed union. A successful downstream recount does not establish the upstream Phasis version.

The following row policy is fixed:

1. Include `core_exact` rows.
2. Include `core_offset` only when no `core_exact` row is present for the same phase, locus, strand, window unit, and expected register.
3. Include `extended_exact` rows.

The physical product identity is `phase + identifier + strand + observed_pos + tag_seq`. Window and register annotations are retained in a many-to-one association ledger; they never multiply a read count. The same tag at different loci remains ambiguous at locus level, while it can be quantified once at global `tag_seq` level if it maps only to Phasis products of one phase.

### Optional hierarchy catalog

The hierarchy metadata is optional. Its TSV schema is:

```text
class    feature_id    parent_level    parent_id
miRNA    miR156a       family          miR156
miRNA    miR156b       family          miR156
```

Without it, a read mapping to `miR156a` and `miR156b` stays in a stable named equivalence group. With it, a separate family-level matrix can count that read once for `miR156`; individual member-level matrices remain conservative.

## Mapping and count assignment

One atomic reference index contains every supplied reference feature and all selected Phasis products. The Bowtie 1 command is equivalent to:

```text
bowtie -f -v 0 -a --best --strata --norc
```

Thus mapping is exact, gap-free, and direct orientation only. Hits at several positions within one supplied feature—for example within a full tRNA—are deduplicated to that one feature. Hits at multiple atomic features become a named equivalence group. Counts are not copied to every hit and are not fractionally allocated in v1.

Primary outputs include:

- global atomic feature counts, including named equivalence groups;
- Phasis 21/24 `tag_seq`, physical-product, locus-plus-tag, and conservative locus-total matrices;
- optional hierarchy-parent matrices;
- per-sequence assignment ledgers, equivalence-group membership, and conservation tables;
- checksums, commands, software versions, catalog decisions, and resolved target metadata.

Window-specific DA is not claimed from sequence-only FASTQ libraries: identical tags with more than one window/register association cannot be distinguished without additional coordinate-level evidence.

## Commands

Validate every file and catalog decision before using Bowtie:

```bash
phasis-srna-da validate \
  --processed-libraries /data/processed_libraries \
  --library-format fastq --min-length 18 --max-length 30 \
  --targets examples/targets.example.tsv \
  --assembly-id Zm-B73-REFERENCE-NAM-5.0 \
  --annotation-id my_frozen_sRNA_reference_v1 \
  --reference-orientation molecule \
  --reference miRNA=/refs/maize_mature_miRNA.fa \
  --reference tRNA=/refs/maize_full_tRNA.fa \
  --reference snRNA=/refs/maize_snRNA.fa \
  --reference snoRNA=/refs/maize_snoRNA.fa \
  --reference rRNA=/refs/maize_rRNA.fa \
  --phasi-catalog /data/21_phasiRNAs.tsv \
  --phasi-catalog /data/24_phasiRNAs.tsv
```

Run mapping and create matrices in a new empty output directory:

```bash
phasis-srna-da run \
  --processed-libraries /data/processed_libraries \
  --library-format fastq --min-length 18 --max-length 30 \
  --targets examples/targets.example.tsv \
  --assembly-id Zm-B73-REFERENCE-NAM-5.0 \
  --annotation-id my_frozen_sRNA_reference_v1 \
  --reference-orientation molecule \
  --reference miRNA=/refs/maize_mature_miRNA.fa \
  --reference tRNA=/refs/maize_full_tRNA.fa \
  --reference snRNA=/refs/maize_snRNA.fa \
  --reference snoRNA=/refs/maize_snoRNA.fa \
  --reference rRNA=/refs/maize_rRNA.fa \
  --phasi-catalog /data/21_phasiRNAs.tsv \
  --phasi-catalog /data/24_phasiRNAs.tsv \
  --hierarchy examples/hierarchy.example.tsv \
  --threads 16 \
  --outdir results/2026-07-24_counting
```

The paths above are placeholders. Substitute existing inputs and use the same reference arguments for validation and counting. On Hive, execute these commands in a compute job and set threads from `SLURM_CPUS_PER_TASK`.

Run DA from the completed count directory without remapping:

```bash
phasis-srna-da da \
  --count-run results/2026-07-24_counting \
  --normalization total-qc-depth \
  --contrast Lobe,Epi \
  --min-count 10 \
  --prefilter-policy all-replicates-in-one-condition \
  --fdr 0.05 --threads 4 \
  --outdir results/2026-07-24_da
```

Condition names must exactly match the target sheet. Repeat `--contrast NUMERATOR,DENOMINATOR` for every intended comparison; positive log2 fold change means higher relative abundance in the numerator. No automatic undeclared pairwise comparisons are run. Alternatively, add `--run-da`, `--normalization total-qc-depth`, and explicit contrasts to `run`.

DA uses raw integer counts and PyDESeq2 negative-binomial models with one fixed size factor per biological sample: retained library abundance divided by the geometric mean across samples. The same factors apply to every projection, including unassigned and rRNA abundance in the depth denominator. This is the currently supported CLI normalization policy; it is not a claim that library composition is biologically invariant. CP30M values are descriptive outputs, not model inputs.

The default prefilter requires at least `--min-count` in every biological replicate of at least one condition. `--prefilter-policy any-samples` is an explicit alternative. PyDESeq2 independent filtering and Cook's outlier handling can leave adjusted p-values missing. Benjamini–Hochberg adjustment is within each projection and contrast. Inspect `feature_filtering.tsv`, `model_metadata.json`, `normalization_factors.tsv`, and any `SKIPPED.txt`; command completion alone does not prove that every model was fitted. A parametric dispersion trend can fall back to a mean trend, recorded in model metadata.

The five default projections are `global_atomic`, `phasi_tag_seq`, `phasi_product`, `phasi_locus_tag`, and `phasi_locus`. Each is a different resolution of the same reads. A locus summary aggregates declared products; it does not measure intact precursor RNA or all reads in a PHAS interval. The CLI does not automatically generate publication figures or a biological interpretation report.
