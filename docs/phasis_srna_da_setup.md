# Set up phasis-srna-da on Hive

This guide prepares the software used in the [small-RNA differential-abundance tutorial](ms28_cli_tutorial.md). `phasis-srna-da` is a reusable downstream command-line tool. You supply processed libraries, a sample sheet, reference FASTAs, and mandatory Phasis product files. The ms28 dataset is an optional worked example.

Installation and analysis run on a **compute node**. If you already have an `srun` shell, start at section 2. All commands here are for you to execute; installation does not start an analysis.

## 1. Optional: connect to Hive and open an interactive allocation

From your own computer, replace `YOUR_UCD_USERNAME` with your campus account name:

```bash
ssh YOUR_UCD_USERNAME@hive.hpc.ucdavis.edu
```

On the login node, check which accounts you can use and start a `screen` session named `phasis-srna-da` so your session can survive a disconnected terminal. Note the login hostname so you can reconnect to the same node later:

```bash
source /etc/profile.d/modules.sh
sacctmgr show assoc user="$USER" format=account%20,partition%20,qos%40
hostname
screen -S phasis-srna-da
```

Inside `screen`, choose an account from your association listing and request a compute shell. This example uses the Genome Center account; replace it if your account access differs:

```bash
SRNA_ACCOUNT=genome-center-grp
srun --account="$SRNA_ACCOUNT" --partition=high \
  --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=64G \
  --time=08:00:00 --pty /bin/bash -l
```

Wait for the allocation. The command requests one task on one node, eight CPUs, 64 GB of memory, and an eight-hour limit. That is an example allocation, not a measured runtime or memory guarantee for your dataset. Set resources according to your workload and account limits. If you use `bcmeyersgrp`, replace the account value. An explicit `--qos` is usually unnecessary: use only a QOS authorized for your selected account/partition.

`high` is non-preemptible. To detach from `screen`, press **Ctrl+A**, then **D**. After reconnecting to the same Hive login node, list sessions with `screen -ls` and resume with `screen -r phasis-srna-da`. Detaching does not extend your allocation's time limit.

## 2. Confirm your compute shell and set your workspace

In the shell opened by `srun`:

```bash
hostname
printf 'Job: %s\nCPUs per task: %s\n' "${SLURM_JOB_ID:-unset}" "${SLURM_CPUS_PER_TASK:-1}"
```

Confirm that the hostname is a compute node, not `login1` or `login2`, and that a Slurm job ID is present.

Choose a group-storage location you can write to. Meyers lab users can use:

```bash
export SRNA_HOME="/quobyte/bcmeyersgrp/$USER/srna-da"
export SRNA_CODE="$SRNA_HOME/software/phasis-srna-da"
export SRNA_ENV="$SRNA_HOME/envs/phasis-srna-da"
export PIP_CACHE_DIR="$SRNA_HOME/cache/pip"
mkdir -p "$SRNA_HOME/software" "$SRNA_HOME/envs" "$PIP_CACHE_DIR" "$SRNA_HOME/studies"
```

`$USER` expands to your Hive username. Change the group-storage prefix if you belong to another lab. Your storage group and Slurm accounting group need not be the same. Environments, caches, and analysis data belong on group storage rather than your 20-GB home directory.

## 3. Download the software from GitHub

The standalone software repository is [thalescherubino/phasis-srna-da](https://github.com/thalescherubino/phasis-srna-da), extracted from the Maize Anther Primordia Spatial Small-RNA Analysis project. Use the development release `v0.2.0.dev0` for this tutorial. Public HTTPS downloads require no GitHub account or token.

Choose an unused `SRNA_CODE` directory, then run:

```bash
git clone --branch v0.2.0.dev0 --single-branch \
  https://github.com/thalescherubino/phasis-srna-da.git "$SRNA_CODE"
git -C "$SRNA_CODE" rev-parse HEAD
ls "$SRNA_CODE/pyproject.toml" "$SRNA_CODE/src/phasis_srna_da/cli.py" \
  "$SRNA_CODE/examples/srna_da_config.example.sh" \
  "$SRNA_CODE/examples/ms28_da_config.example.sh"
```

Continue only if the clone succeeds. Checking out a release tag may show a normal “detached HEAD” message. Keep this checkout unchanged during the analysis and record its commit below. If you already have this release, set `SRNA_CODE` to its location and skip cloning.

Installing the software does not download the annotation or example libraries. Obtain them separately using the [data-download guide](phasis_srna_da_downloads.md). This is a development release; software tests do not establish successful completion or biological validity of the full ms28 tutorial.

## 4. Load Python and Bowtie 1

Still on the compute node:

```bash
source /etc/profile.d/modules.sh
module avail python/3.11.9 bowtie/1.3.1
module load python/3.11.9 bowtie/1.3.1
python --version
bowtie --version
bowtie-build --version
```

Python must be at least 3.10. These instructions use Hive's Python 3.11.9 and Bowtie 1.3.1 modules. **Bowtie 1 is required**; Bowtie 2 is not a substitute. If module names change, inspect `module avail` / `module search` before choosing replacements.

## 5. Create a separate Python environment and install the tool

For a first installation, choose an unused `SRNA_ENV` path, then run:

```bash
python -m venv "$SRNA_ENV"
export SRNA_PYTHON="$SRNA_ENV/bin/python"
"$SRNA_PYTHON" -m pip install --upgrade pip
"$SRNA_PYTHON" -m pip install "${SRNA_CODE}[da]"
```

`[da]` installs the optional DA dependencies, including NumPy and PyDESeq2, in addition to the base package dependencies. The package currently requires PyDESeq2 `>=0.5,<0.6`. Pip uses its normal isolated build environment to provide setuptools. Keep build isolation enabled.

This installs a local copy of the software; it is not a Phasis installation and does not change the environment used for upstream discovery. The commands use the full path to the environment's Python, so activation is optional. Subsequent source edits require reinstalling if you intend to use them; keep one installed version unchanged throughout an analysis.

If you already have a separate environment for this tool, skip environment creation and set `SRNA_ENV` and `SRNA_PYTHON` to it before verification.

## 6. Verify the installation and save its identity

```bash
"$SRNA_PYTHON" -m pip check
"$SRNA_PYTHON" -m phasis_srna_da --version
"$SRNA_PYTHON" -m phasis_srna_da --help
"$SRNA_PYTHON" - <<'PY'
import sys
from importlib.metadata import version
import numpy
import pandas
from pydeseq2.dds import DeseqDataSet
import phasis_srna_da

print('Python:', sys.version.split()[0])
print('Loaded module:', phasis_srna_da.__file__)
for name in ('phasis-srna-da', 'pydeseq2', 'numpy', 'pandas'):
    print(name + ':', version(name))
assert phasis_srna_da.__version__ == version('phasis-srna-da'), 'Source/package version mismatch'
PY
```

Expect `pip check` to report no broken requirements, CLI help to list `validate`, `run`, and `da`, and the imports/version comparison to succeed. Fix installation errors before running validation on your libraries.

Record the installation on group storage:

```bash
SRNA_SETUP_RECORD="$SRNA_HOME/software/setup_$(date +%Y%m%d_%H%M%S)"
mkdir "$SRNA_SETUP_RECORD"
git -C "$SRNA_CODE" rev-parse HEAD > "$SRNA_SETUP_RECORD/source_commit.txt"
"$SRNA_PYTHON" -m pip freeze > "$SRNA_SETUP_RECORD/python_packages.txt"
"$SRNA_PYTHON" -m phasis_srna_da --version > "$SRNA_SETUP_RECORD/cli_version.txt"
bowtie --version > "$SRNA_SETUP_RECORD/bowtie_version.txt"
sha256sum "$SRNA_CODE/pyproject.toml" "$SRNA_CODE"/src/phasis_srna_da/*.py \
  > "$SRNA_SETUP_RECORD/source.sha256"
```

These records identify your installation; they are not a fully portable dependency lockfile. Keep the source copy and use the recorded versions when preparing a reproducible lab release.

## 7. Reuse the environment in a later session

In a new `srun` compute shell, restore the paths and load the same modules; there is no need to reinstall:

```bash
export SRNA_HOME="/quobyte/bcmeyersgrp/$USER/srna-da"
export SRNA_CODE="$SRNA_HOME/software/phasis-srna-da"
export SRNA_ENV="$SRNA_HOME/envs/phasis-srna-da"
export SRNA_PYTHON="$SRNA_ENV/bin/python"
export PIP_CACHE_DIR="$SRNA_HOME/cache/pip"
source /etc/profile.d/modules.sh
module load python/3.11.9 bowtie/1.3.1
"$SRNA_PYTHON" -m phasis_srna_da --version
```

Adjust those paths if you chose different locations. Obtain the [B73v5 annotation and example inputs](phasis_srna_da_downloads.md), then continue with the [analysis tutorial](ms28_cli_tutorial.md), choosing either your own dataset or the ms28 example. For cluster account or administrative software problems, contact `hpc-help@ucdavis.edu`.
