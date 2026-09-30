# Small-RNA differential abundance with phasis-srna-da

`phasis-srna-da` is a reusable command-line workflow for counting processed small-RNA libraries against a fixed reference and testing relative differential abundance. This tutorial explains how to use your own dataset and provides **B73 WT/ms28 anthers as a worked example**. The ms28 example is distinct from the anther-primordia spatial dataset.

**Start with [software setup on Hive](phasis_srna_da_setup.md).** It installs the tool in one Conda environment and loads Bowtie 1.

The analysis below starts **inside an existing `srun` compute-node shell**. You run each stage yourself: **prepare inputs → validate → count → review QC → DA → inspect results**. No step submits another job automatically.

## 1. What you need for any dataset

| Input or decision | Required contents |
|---|---|
| Processed libraries | One `.tag` or `.fastq.gz` file per sample/library, using one declared format per run |
| Target sheet | Tab-separated sample names, condition labels, and biological replicate identifiers |
| Class reference FASTAs | At least one explicitly labelled class; include all classes you intend to quantify |
| Phasis product files | One or more fixed `*_phasiRNAs.tsv` catalogs; mandatory even with custom non-PHAS references |
| Reference identities | Assembly ID and stable frozen annotation/reference release ID |
| Read lengths | Explicit minimum and maximum retained lengths |
| Comparisons | Explicit numerator and denominator condition names matching the target sheet |

Use the same upstream processing across samples. The tool does not trim adapters, perform quality trimming, discover loci, or retrieve libraries from a database. A class removed from the input libraries cannot be recovered by adding its reference later.

### Library formats

For `--library-format tag`, use a two-column file without a header:

```text
ACGTACGTACGTACGTACGTA	120
TGCATGCATGCATGCATGCAT	35
```

The separator must be an actual tab. Each row contributes its positive integer abundance. Identical canonical sequences are summed.

For `--library-format fastq`, supply gzip-compressed, four-line FASTQ records, such as:

```text
@read_1
ACGTACGTACGTACGTACGTA
+
IIIIIIIIIIIIIIIIIIIII
```

Each retained FASTQ record contributes one count; abundance is not inferred from its header. Complete records containing **N** are discarded before length filtering. Other ambiguous bases are rejected. Case is normalized and U is converted to T. These rules apply to custom datasets as well as the example.

### Target sheet

A custom study with two conditions and three replicates could use:

```text
sample_id	condition	biological_replicate	library_id
control_r1	control	1	control_1
control_r2	control	2	control_2
control_r3	control	3	control_3
treatment_r1	treatment	1	treatment_1
treatment_r2	treatment	2	treatment_2
treatment_r3	treatment	3	treatment_3
```

Save this as a tab-separated `targets.tsv`. Here `library_id=control_1` resolves to `control_1.tag` or `control_1.fastq.gz`, according to your chosen format. If filenames already match `sample_id`, omit `library_id`. It is a filename stem, not a requirement to use external database identifiers.

Use at least two independent biological replicates per condition for inferential DA; three are recommended. Technical runs are not additional biological replicates. Optional `technical_replicate` labels allow explicitly declared technical libraries from the same biological replicate to be combined. `batch` is optional; incomplete or confounded batch information cannot provide a valid batch adjustment.

### Reference and Phasis files

Provide one `--reference CLASS=FASTA` option for each class. FASTA feature IDs must be unique within a class. References must contain unambiguous sequences in **5′-to-3′ molecular orientation**. `--reference-orientation molecule` declares this contract; it does not convert genomic-strand exports automatically.

Mandatory Phasis TSV columns are `identifier`, `phase`, `window_unit_id`, `strand`, `expected_register_pos`, `register_class`, and `tag_seq`. Preserve the original export's additional fields, including observed positions and provenance. Supported phases are 21 and 24, with tag length matching phase. Do not manufacture Phasis product files from intervals or generic FASTAs.

The catalog includes `core_exact` and `extended_exact` rows. A `core_offset` row is included only when its corresponding exact core register is absent. Discovery abundances in the catalog do not become sample counts: the analysis recounts your supplied libraries. Use a single fixed catalog per phase unless you have an explicitly documented curated union. Your inputs must share compatible assembly and annotation identities; supplying an ID is a declaration, not an automatic compatibility test.

## 2. Restore your software paths

If you continued directly from the [setup guide](phasis_srna_da_setup.md), these paths are already defined. Otherwise restore them inside your compute shell:

```bash
export SRNA_HOME="/quobyte/bcmeyersgrp/$USER/srna-da"
export SRNA_CODE="$SRNA_HOME/software/phasis-srna-da"
export SRNA_ENV="$SRNA_HOME/conda/phasis-srna-da"
source /etc/profile.d/modules.sh
module load conda bowtie/1.3.1
conda activate "$SRNA_ENV"
hostname
printf 'Job: %s\nThreads: %s\n' "${SLURM_JOB_ID:-unset}" "${SLURM_CPUS_PER_TASK:-1}"
phasis-srna-da --version
```

Change the group-storage prefix if needed. `$USER` is appropriate for **your own workspace**. Shared libraries and frozen references may belong to another owner: use their actual supplied paths rather than replacing the owner's name with `$USER`.

Use the same environment throughout an analysis. The commands below take their thread count from `SLURM_CPUS_PER_TASK`, falling back to one thread. All validation, counting, model fitting, and bulk checksum work stays in the compute session.

## 3. Choose your dataset configuration

Choose **one** of the following routes. Both lead to the same validation, counting, and DA commands.

### A. Your own dataset

Create a study directory and copy the generic template:

```bash
export SRNA_STUDY="$SRNA_HOME/studies/my_study"
mkdir -p "$SRNA_STUDY"
cp "$SRNA_CODE/examples/srna_da_config.example.sh" "$SRNA_STUDY/analysis_config.sh"
```

Open `analysis_config.sh` in your editor and replace every placeholder. The [generic template](../examples/srna_da_config.example.sh) defines two Bash arrays: `SRNA_INPUT_ARGS` for validation/counting and `SRNA_DA_ARGS` for the model settings and contrasts. It only defines arguments; sourcing it does not run analysis.

Edit the following:

| Setting | Your replacement |
|---|---|
| `SRNA_LIBRARIES` | Directory containing your processed files |
| `SRNA_TARGETS` | Your target sheet |
| `SRNA_REFERENCES`, `SRNA_PHASIS` | Your fixed FASTA and Phasis product locations |
| `--library-format` | `tag` or `fastq` |
| `--min-length`, `--max-length` | Length bounds justified for your study |
| `--assembly-id`, `--annotation-id` | Actual reference identities; replace the `REPLACE_WITH_...` labels |
| Repeated `--reference` entries | Exact class labels and FASTA filenames you supply |
| Repeated `--phasiRNAs` entries | Your mandatory Phasis file(s); remove a phase only if it is not part of your fixed catalog |
| Repeated `--contrast` entries | Actual condition labels, written `NUMERATOR,DENOMINATOR` |

Remove FASTA entries for classes you do not supply and add entries for other declared classes. Keep at least one reference FASTA and one Phasis catalog. The generic template's 18–50-nt bounds and `treatment,control` contrast are examples to review, not assumptions about your study. Optional hierarchy metadata can be supplied with `--hierarchy /path/to/hierarchy.tsv`; see the [input contract](phasis_srna_da.md).

### B. The ms28 worked example

This example uses 18 processed B73 WT/ms28 anther libraries: two genotypes, three stages (0.4, 2, and 5 mm), and three declared replicates per group. They are **workflow-test data, not reliable biological validation**.

```bash
export SRNA_STUDY="$SRNA_HOME/studies/ms28_example"
mkdir -p "$SRNA_STUDY"
cp "$SRNA_CODE/examples/ms28_da_config.example.sh" "$SRNA_STUDY/analysis_config.sh"
```

Edit `SRNA_EXAMPLE_ROOT` in the copied file to your example location. The layout below is a portable template; if your existing inputs use other paths, edit the configuration to those paths without moving or changing the frozen originals. The [ms28 template](../examples/ms28_da_config.example.sh) expects these paths relative to that location:

| Input | Relative location |
|---|---|
| Tag libraries | `libraries/` |
| Target sheet | `targets.tsv` |
| Frozen FASTAs | `annotation/fasta/` |
| Mandatory Phasis products | `phasiRNAs/` |
| Exploratory PHAS-like FASTAs | `candidate_fasta/` |

These data are supplied separately from the software. Follow the [data-download guide](phasis_srna_da_downloads.md) for the public B73v5 reference bundle. The processed libraries and target sheet are supplied separately and must be verified before use. Other result/report packages do not replace the tutorial input bundle.

The example root may be a shared directory; the analysis reads it without rewriting the reference or libraries. If the lab distributes a differently organized bundle, edit the individual paths in your configuration. Preserve the reference and library inputs separately.

The template retains 18–50 nt, discards N-containing records, and includes miRNA, tRNA, rRNA, snRNA, snoRNA, 21-/24-nt phasiRNA, unsupported MIR candidates, and 21-/24-PHAS-like candidates. Candidate labels remain exploratory. The frozen assembly is `Zm-B73-REFERENCE-NAM-5.0`; the annotation identity is `B73v5_curated_sRNA_reference_v1_2026-08-27__PHAS_molecule_oriented_v1_2026-08-30`.

The nine example comparisons are:

| Family | Numerator versus denominator |
|---|---|
| Genotype within stage | `ms28_0p4,WT_0p4`; `ms28_2,WT_2`; `ms28_5,WT_5` |
| Stage within WT | `WT_2,WT_0p4`; `WT_5,WT_0p4`; `WT_5,WT_2` |
| Stage within ms28 | `ms28_2,ms28_0p4`; `ms28_5,ms28_0p4`; `ms28_5,ms28_2` |

These are pairwise condition comparisons, not a genotype-by-stage interaction test.

## 4. Start a new run and save the settings

After editing your chosen configuration:

```bash
bash -n "$SRNA_STUDY/analysis_config.sh"
cat "$SRNA_STUDY/analysis_config.sh"

SRNA_RUN="$SRNA_STUDY/results/run_01"
mkdir -p "$SRNA_STUDY/results"
mkdir "$SRNA_RUN"
mkdir "$SRNA_RUN/logs"
cp "$SRNA_STUDY/analysis_config.sh" "$SRNA_RUN/analysis_config.sh"
source "$SRNA_RUN/analysis_config.sh"
set -o pipefail

python -m pip freeze > "$SRNA_RUN/python_packages.txt"
phasis-srna-da --version > "$SRNA_RUN/cli_version.txt"
bowtie --version > "$SRNA_RUN/bowtie_version.txt"
sha256sum "$SRNA_RUN/analysis_config.sh" "$SRNA_CODE/pyproject.toml" \
  "$SRNA_CODE"/src/phasis_srna_da/*.py > "$SRNA_RUN/software_and_config.sha256"
```

Choose an unused run name. If `mkdir "$SRNA_RUN"` reports that it exists, stop and choose `run_02`, etc.; do not overwrite a prior run's logs. Keep the saved configuration and installed software unchanged through all stages. `pipefail` ensures a failed CLI command is not hidden by a successful `tee` log-writing process.

## 5. Validate the inputs

```bash
phasis-srna-da validate \
  "${SRNA_INPUT_ARGS[@]}" --json \
  > "$SRNA_RUN/validation.json" 2> "$SRNA_RUN/logs/validation.err"
```

Wait for the command to return. If it reports an error, inspect `logs/validation.err` and fix the inputs before continuing. On success:

```bash
cat "$SRNA_RUN/validation.json"
```

Check `valid: true`, the sample/condition counts, reference classes, Phasis file count, assembly/annotation identities, length bounds, orientation, N exclusions and retained abundance. Expected counts come from **your** design. For ms28, expect 18 samples, six conditions, two Phasis product files, and ten class labels.

Validation reads the complete libraries but does not map them or fit DA models. Review warnings before proceeding.

## 6. Count against the fixed reference

```bash
phasis-srna-da run \
  "${SRNA_INPUT_ARGS[@]}" \
  --threads "${SLURM_CPUS_PER_TASK:-1}" \
  --outdir "$SRNA_RUN/counts_run" \
  2>&1 | tee "$SRNA_RUN/logs/count.log"
```

The foreground command reports library progress. All supplied classes compete in one Bowtie 1 index, using direct perfect matches: `-v 0 -a --best --strata --norc`. A Phasis product requires full-sequence identity. Generic FASTA features allow substrings; their counts represent reference-associated signal rather than a universal exact-mature-product measurement.

Each retained abundance is assigned once within a projection: a unique unit, a named equivalence group, or unassigned. Multiple matches do not multiply counts or receive fractionally rounded allocations.

After the command succeeds:

```bash
cat "$SRNA_RUN/counts_run/provenance/run_metadata.json"
column -t -s $'\t' "$SRNA_RUN/counts_run/qc/library_summary.tsv" | less -S
column -t -s $'\t' "$SRNA_RUN/counts_run/mapping/conservation_by_sample.tsv" | less -S
```

Press `q` to leave `less`. Require `status: complete`, one QC row per declared library, and `conserved=true` for every sample/projection:

```text
unique assigned + equivalence-group assigned + unassigned = retained abundance
```

| Matrix under `counts_run/counts/` | Measurement |
|---|---|
| `global_atomic.raw.tsv` | All supplied features and equivalence groups |
| `phasi_tag_seq.raw.tsv` | PhasiRNA sequence units |
| `phasi_product.raw.tsv` | Declared physical Phasis products |
| `phasi_locus_tag.raw.tsv` | Locus-plus-sequence units |
| `phasi_locus.raw.tsv` | Conservative locus output from declared products |

These are different resolutions of the same reads. A locus summary does not measure intact precursor RNA or automatically count every read inside a PHAS interval. Review counting/QC before starting DA.

## 7. Run the declared differential-abundance comparisons

```bash
phasis-srna-da da \
  --count-run "$SRNA_RUN/counts_run" \
  "${SRNA_DA_ARGS[@]}" \
  --threads "${SLURM_CPUS_PER_TASK:-1}" \
  --outdir "$SRNA_RUN/da_run" \
  2>&1 | tee "$SRNA_RUN/logs/da.log"
```

This uses completed integer counts without remapping. It checks count-run evidence and fits each projection separately. Positive log2 fold change means higher relative abundance in the numerator condition.

The supplied templates declare `total-qc-depth`: each sample's retained abundance divided by the geometric mean depth gives its size factor, shared across projections. Retained rRNA and unassigned abundance are included in this denominator. The default prefilter requires at least 10 raw counts in every biological replicate of at least one condition. FDR is 0.05 within each projection and contrast. These are explicit settings to review for your study, not evidence that composition or total RNA abundance is biologically invariant.

After the command succeeds:

```bash
cat "$SRNA_RUN/da_run/provenance/run_metadata.json"
ls "$SRNA_RUN/da_run/da/phasi_tag_seq/"
cat "$SRNA_RUN/da_run/da/phasi_tag_seq/model_metadata.json"
```

Inspect every projection under `da_run/da/`:

| Output | Review |
|---|---|
| `feature_filtering.tsv` | Units retained/rejected by the raw-count prefilter |
| `normalization_factors.tsv` | Sample depths and size factors |
| `model_metadata.json` | Design, contrasts, dispersion fit and fallback |
| `NUMERATOR_vs_DENOMINATOR.tsv` | Fold changes, p-values, adjusted p-values, feature labels |
| `SKIPPED.txt`, if present | Reason that model was withheld |

For a custom study, expect one result table per declared contrast **for each fitted projection**. The ms28 example produces 45 tables only if all five default projections fit. A successful process exit alone does not establish that every model ran. Independent filtering and Cook's outlier handling can leave missing p-values/adjusted p-values. Record skipped models and dispersion fallbacks. CP30M tables are descriptive; fitting uses raw integer counts.

## 8. Preserve evidence and confirm input integrity

Before treating a result as current, check that source inputs still match the hashes recorded during counting. Run this in the compute session:

```bash
python - "$SRNA_RUN/counts_run/provenance/input_manifest.tsv" <<'PY'
import csv
from pathlib import Path
import sys
from phasis_srna_da.io import sha256_file

with open(sys.argv[1]) as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        if sha256_file(Path(row['path'])) != row['sha256']:
            raise SystemExit('Input changed: ' + row['path'])
print('All count-run input checksums still match.')
PY
```

Retain the complete run and logs. Copy small provenance records to your project's tracked manifest area, for example:

```bash
SRNA_RECORD="$SRNA_STUDY/provenance/$(basename "$SRNA_RUN")"
mkdir -p "$SRNA_STUDY/provenance"
mkdir "$SRNA_RECORD"
cp "$SRNA_RUN/analysis_config.sh" "$SRNA_RUN/python_packages.txt" \
  "$SRNA_RUN/cli_version.txt" "$SRNA_RUN/bowtie_version.txt" \
  "$SRNA_RUN/software_and_config.sha256" "$SRNA_RECORD/"
cp "$SRNA_RUN/counts_run/provenance/input_manifest.tsv" "$SRNA_RECORD/count_inputs.tsv"
cp "$SRNA_RUN/counts_run/provenance/run_metadata.json" "$SRNA_RECORD/count_metadata.json"
cp "$SRNA_RUN/counts_run/provenance/commands.tsv" "$SRNA_RECORD/mapping_commands.tsv"
cp "$SRNA_RUN/da_run/provenance/input_manifest.tsv" "$SRNA_RECORD/da_inputs.tsv"
cp "$SRNA_RUN/da_run/provenance/run_metadata.json" "$SRNA_RECORD/da_metadata.json"
for SRNA_RESULT in "$SRNA_RUN/counts_run/counts/"*.raw.tsv \
  "$SRNA_RUN/da_run/da/"*/*_vs_*.tsv; do
  if [[ -f "$SRNA_RESULT" ]]; then sha256sum "$SRNA_RESULT"; fi
done > "$SRNA_RECORD/result_checksums.sha256"
```

Include those small records in your normal version-control review. They retain source paths, checksums, assembly/annotation identities, commands, versions, parameters, and status. Keep libraries, reference sequences, indexes, alignments, bulk matrices, plots, environments and caches out of Git. A provenance directory is not automatically tracked merely because it exists.

## 9. Stop, reconnect, or troubleshoot

Press **Ctrl+C** to interrupt the foreground command. Preserve incomplete outputs and logs; a directory's existence does not mean a stage completed. Counting cannot resume midway through a library: use a new run/output directory for an incomplete count attempt. Completed counts can be reused for DA after restoring the same environment and paths; use a new DA output directory if an earlier DA attempt left partial files.

| Situation | Next action |
|---|---|
| Missing module or Python dependency | Return to the setup guide and verify the intended environment |
| Missing library | Check `library_id`/`sample_id`, basename, suffix, and `--library-format` |
| Validation or mapping error | Inspect the saved error/log file before running later stages |
| Input checksum mismatch | Investigate the changed source and preserve the frozen reference |
| Existing nonempty output directory | Choose a new output directory |
| DA `SKIPPED.txt` | Read the design/feature-count reason |
| Allocation ends or exceeds memory | Obtain suitable resources separately; preserve logs before retrying |

Inspect an allocation with `squeue -j JOBID` or `sacct -j JOBID --format=JobID,State,ExitCode,Elapsed,MaxRSS`. `scancel JOBID` releases the entire allocation and stops work in it. Run one writer per output directory. Administrative problems go to `hpc-help@ucdavis.edu`.

The output is a downstream recount, not a new Phasis discovery or proof of an upstream version. Reads in a PHAS interval are not automatically phasiRNAs, a 24-nt-rich interval is not automatically 24-PHAS, and spatial precursor/mature patterns alone do not establish transport. For the ms28 example, assess computational behavior rather than treating results as reliable biological validation.
