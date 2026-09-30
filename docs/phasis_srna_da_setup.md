# Set up `phasis-srna-da` on Hive

Use these commands inside an existing `srun` compute shell. They create one
Conda environment and install the public command-line tool. They do not start
an analysis.

## Install once

```bash
source /etc/profile.d/modules.sh
module load conda bowtie/1.3.1

export SRNA_HOME="/quobyte/bcmeyersgrp/$USER/srna-da"
export SRNA_ENV="$SRNA_HOME/conda/phasis-srna-da"
export SRNA_CODE="$SRNA_HOME/software/phasis-srna-da"
export CONDA_PKGS_DIRS="$SRNA_HOME/conda-pkgs"
export PIP_CACHE_DIR="$SRNA_HOME/pip-cache"
mkdir -p "$SRNA_HOME/conda" "$SRNA_HOME/software" "$CONDA_PKGS_DIRS" "$PIP_CACHE_DIR"

conda create --yes --no-default-packages --prefix "$SRNA_ENV" python=3.11 pip
conda activate "$SRNA_ENV"

git clone --branch v0.2.0.dev0 --single-branch \
  https://github.com/thalescherubino/phasis-srna-da.git "$SRNA_CODE"
python -m pip install "${SRNA_CODE}[da]"

phasis-srna-da --version
bowtie --version
```

Keep `SRNA_HOME` on Quobyte so the environment, package cache, software, and
downloaded inputs are available in later sessions.

## Use it later

Start the `srun` shell as usual, then run:

```bash
source /etc/profile.d/modules.sh
module load conda bowtie/1.3.1

export SRNA_HOME="/quobyte/bcmeyersgrp/$USER/srna-da"
export SRNA_ENV="$SRNA_HOME/conda/phasis-srna-da"
conda activate "$SRNA_ENV"

phasis-srna-da --help
```

Next, follow [the reference-download guide](phasis_srna_da_downloads.md), then
[the CLI tutorial](ms28_cli_tutorial.md).
