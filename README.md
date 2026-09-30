# phasis-srna-da

Conservative counting and relative differential abundance of small RNAs using fixed Phasis product catalogs and class-labelled reference FASTAs.

The workflow accepts processed `.fastq.gz` or two-column `.tag` libraries, discards records containing N, and uses explicitly declared read-length bounds. All supplied classes compete in one Bowtie 1 index. Multi-mapping abundance is retained in named equivalence groups. Fixed Phasis TSV inputs are mandatory; the software does not run Phasis discovery.

This standalone tool originated in the Maize Anther Primordia Spatial Small-RNA Analysis project and supports other datasets. The optional B73 WT/ms28 example is separate from the spatial project and is intended as workflow-test data.

## Download and install

Requirements: Python 3.10 or later, Git, Bowtie **1**, and `bowtie-build`. On Hive, run installation and analysis inside your existing `srun` compute shell and keep the environment on group storage. The [Hive setup guide](docs/phasis_srna_da_setup.md) includes modules, storage paths and an optional named `screen` session.

From a writable software directory:

```bash
git clone --branch v0.2.0.dev0 --single-branch \
  https://github.com/thalescherubino/phasis-srna-da.git
cd phasis-srna-da
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install '.[da]'
.venv/bin/python -m pip check
.venv/bin/python -m phasis_srna_da --version
.venv/bin/python -m phasis_srna_da --help
bowtie --version
bowtie-build --version
```

`[da]` installs PyDESeq2 and the other model dependencies. A normal installation copies the package into the environment; reinstall after changing source. Use one fixed version per analysis. A GitHub account is not required to clone this public repository.

## Prepare and run your analysis

Follow the [step-by-step CLI tutorial](docs/ms28_cli_tutorial.md): prepare inputs, validate, count, review QC, run explicit DA contrasts, and inspect results. The CLI stages are `validate`, `run`, and `da`.

- [Required inputs and counting/model contract](docs/phasis_srna_da.md)
- [Generic study configuration](examples/srna_da_config.example.sh)
- [B73 WT/ms28 example configuration](examples/ms28_da_config.example.sh)
- [B73v5 annotation download and example-data availability](docs/phasis_srna_da_downloads.md)

The software repository does not include sequence libraries, biological reference data, indexes, or bulk results. The public B73v5 annotation alone does not supply every input for the ms28 example. The exact processed libraries still need a public download bundle, and the molecule-oriented Phasis products must be supplied separately.

## Development status and testing

`0.2.0.dev0` is a development release. GitHub Actions installs the package and DA dependencies, checks wheel/source builds, and runs the included tests on synthetic inputs. These checks do not establish biological validity or successful completion of the 18-library ms28 workflow. The existing counting test uses a simulated Bowtie executable; real-data mapping remains a separate validation step.

To run the software tests in a suitable compute environment:

```bash
.venv/bin/python -m pip install '.[da,dev]'
.venv/bin/python -m pytest
```

The initial source extraction is documented in [release notes](RELEASE_NOTES.md) and the [software checksum manifest](provenance/software_manifest.tsv). Reports and dataset-specific cluster jobs from the parent project are not part of the standalone package.

## Interpretation and license

Results measure relative abundance. Reads within a PHAS interval are not automatically phasiRNAs; counts associated with a generic reference sequence are not universally exact mature-product measurements. A downstream recount does not establish the upstream Phasis version.

The software is distributed under the [MIT license](LICENSE). Separately distributed biological data retain their own terms and provenance.
