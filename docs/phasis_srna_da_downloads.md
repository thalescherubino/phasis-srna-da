# Obtain the B73v5 reference and tutorial inputs

Complete [software setup](phasis_srna_da_setup.md) first. The commands here are for you to execute inside your existing `srun` compute shell, with `SRNA_HOME` and `SRNA_PYTHON` defined. They download data and check integrity; they do not run counting or differential abundance.

## 1. Available package and access check

The supplied [Google Drive folder](https://drive.google.com/drive/folders/1uZx6WmrwaZSc1OK9MVkL1LsXanszFG72) contains `B73v5_curated_sRNA_reference_v1_2026-08-27`. Its folder listings and README were accessible without signing in on 2026-09-30. A complete recursive download and full reference checksum verification have **not** been performed as part of preparing this guide.

The reference README identifies assembly `Zm-B73-REFERENCE-NAM-5.0` and describes these components:

| Directory | Contents |
|---|---|
| `annotation/` | Combined and per-class GFF3 annotations |
| `fasta/` | Canonical and separately labelled exploratory reference sequences |
| `tables/` | Feature metadata, equivalence groups, locus and observation tables |
| `classification_changes/` | Curation decisions and supporting diagnostics |
| `documentation/` | Reference documentation |
| `provenance/` | Release manifest, validation records and `SHA256SUMS` |

The same Drive root contains additional report material. The commands below select the reference directories, and do not treat other folders or existing DA result tables as tutorial inputs.

## 2. Download and verify the annotation

Use [gdown](https://github.com/wkentaro/gdown), which supports public Drive files and folders. Install it into a separate download environment so the analysis environment stays fixed. Choose an unused environment path for the first installation:

```bash
SRNA_DOWNLOAD_ENV="$SRNA_HOME/envs/gdown"
SRNA_DOWNLOAD_PYTHON="$SRNA_DOWNLOAD_ENV/bin/python"
"$SRNA_PYTHON" -m venv "$SRNA_DOWNLOAD_ENV"
"$SRNA_DOWNLOAD_PYTHON" -m pip install gdown
"$SRNA_DOWNLOAD_PYTHON" -m gdown --version
```

If that environment already exists, reuse its Python and check the version instead of recreating it. Save downloader versions and choose a new download destination outside the software repository:

```bash
SRNA_B73_REFERENCE="$SRNA_HOME/data/B73v5_curated_sRNA_reference_v1_2026-08-27"
mkdir -p "$SRNA_HOME/data" "$SRNA_HOME/download_records"
"$SRNA_DOWNLOAD_PYTHON" -m pip freeze \
  > "$SRNA_HOME/download_records/gdown_packages_$(date +%Y%m%d_%H%M%S).txt"
```

The following block downloads the six reference directories and README by their public IDs. It leaves separate report packages on Drive. `mkdir` requires a new destination; if it exists, choose a fresh path before retrying. The subshell stops on a download or checksum failure.

```bash
(
  set -e
  mkdir "$SRNA_B73_REFERENCE"
  "$SRNA_DOWNLOAD_PYTHON" -m gdown --no-cookies \
    'https://drive.google.com/uc?id=11nWyjUgdkJcq6MkMUhCUuRj8KdXnXxCX' \
    -O "$SRNA_B73_REFERENCE/README.md"
  while read -r SRNA_FOLDER_NAME SRNA_FOLDER_ID; do
    "$SRNA_DOWNLOAD_PYTHON" -m gdown --folder --no-cookies \
      "https://drive.google.com/drive/folders/$SRNA_FOLDER_ID" \
      -O "$SRNA_B73_REFERENCE/$SRNA_FOLDER_NAME"
  done <<'FOLDERS'
annotation 1yGbTXl9pAP6d1RZ6GtOoA7ypCITd16eT
classification_changes 1Lyr20UY85BnbOXRrzWMgSvNEZX1Kyn2x
documentation 1BILw_lC2JAM4jwxth7I-zTs72_SAd6Vt
fasta 1C8VCreH5Xb-FZvwv6Yefin7K_dqKQ238
provenance 1iJixYYy6CTjDnTbfnXR5EaPbdRbTJWFT
tables 12F9OsOmhL50Vzl7IgatkOTLUmQIST9aV
FOLDERS
  cd "$SRNA_B73_REFERENCE"
  sha256sum --check provenance/SHA256SUMS
)
```

Require successful downloads and all listed checksums to report `OK`. Missing files or checksum mismatches must be resolved before analysis. Do not overwrite the frozen source or regenerate its checksum file to make a mismatch pass. The checksum list covers its declared files; additional files on Drive are not automatically members of the frozen release.

If Drive denies access, throttles downloads, or cannot enumerate a folder completely, stop and retain the error. Do not accept a partial folder as a complete reference. You can also download the reference directories in a browser and transfer them to group storage, preserving the directory structure, then run the same checksum check in your compute shell. A future Zenodo mirror should provide a versioned archive and an externally recorded SHA-256 for that archive.

## 3. Connect the download to the analysis

In your copied `analysis_config.sh`, the downloaded canonical FASTAs can be selected with:

```bash
SRNA_REFERENCE="$SRNA_HOME/data/B73v5_curated_sRNA_reference_v1_2026-08-27/fasta"
```

Use `SRNA_REFERENCES` instead for the generic configuration template, which uses that plural variable name. Its example FASTA basenames also need to match the downloaded files: for example, `miRNA_bona_fide.mature_and_star.fa` and `tRNA.CCA_aware_quantification.fa`.

**The annotation package alone is not a complete ms28 CLI example.** Keep these additional inputs distinct:

| Required input | What remains to supply or verify |
|---|---|
| 18 processed ms28/WT libraries and target sheet | Public download bundle for these exact libraries is pending |
| Fixed 21- and 24-nt Phasis TSV products | The template requires `21_PHAS.molecule_oriented.phasiRNAs.tsv` and `24_PHAS.molecule_oriented.phasiRNAs.tsv`; these exact files were not identified in the inspected reference listings |
| Molecule-oriented PHAS-like candidate FASTAs | The template uses the separately frozen `*.molecule_oriented.fa` products; equivalence to the older reference exports has not been established |

The reference's `21_PHAS.phasiRNAs.fa` and `24_PHAS.phasiRNAs.fa` do not replace mandatory Phasis TSVs. Its compressed `*_phasiRNA_observations.tsv.gz` tables are not verified substitutes either. Do not rename, decompress, or relabel them into the required inputs without checking their schema, register-selection policy, molecular orientation, and reference identity. Preserve the frozen annotation unchanged.

Once the additional inputs are available, follow the [analysis tutorial](ms28_cli_tutorial.md). Retain the agreed 18–50-nt range, N exclusion, all supplied classes, and the interpretation of ms28 as workflow-test data.

## 4. Software release and data distribution

The software repository is [thalescherubino/phasis-srna-da](https://github.com/thalescherubino/phasis-srna-da); follow the [setup guide](phasis_srna_da_setup.md) to clone and install `v0.2.0.dev0`. The repository contains the CLI, documentation, configuration templates and software tests. It does not bundle biological inputs.

| Material | Recommended home |
|---|---|
| CLI source, installation metadata, software tests, tutorials, configuration templates | Dedicated public GitHub repository |
| Frozen B73v5 annotation and its provenance/checksums | Versioned Zenodo dataset record; Drive can remain a mirror |
| Complete tutorial input bundle | Separate versioned dataset archive once the exact inputs and public distribution scope are settled |
| Bulk results, indexes, alignments and private libraries | Outside software Git history |

Zenodo supports version-specific records and DOIs, making an explicitly pinned dataset version suitable for the tutorial. Its standard record allowance is 50 GB and 100 files; package a multi-file reference as an archive with a README and checksums. [Versioning documentation](https://help.zenodo.org/docs/deposit/manage-versions/), [file limits](https://help.zenodo.org/docs/deposit/manage-files/).

GitHub Releases can serve an archive as an alternative or mirror; its release documentation specifies an individual asset limit below 2 GiB. Ordinary Git rejects files larger than 100 MiB. Regardless of size, keep reference sequences and bulk data out of this project's Git history. [Release assets](https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases), [Git file limits](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github).

Prepare the annotation archive from the manifest-defined frozen release, preserving its bytes. Do not archive the entire current Drive root as the annotation release: it also holds report material and may include other additions. Keep archive checksums, per-file checksums, source paths, assembly/annotation identities, release commands and version/status records. Publish any required molecule-oriented Phasis products as separately identified inputs without modifying the original reference release.

Add exact Zenodo record/file URLs and archive checksums after deposition. Until then, the Drive instructions above are the available annotation download route. Full public reproduction of the 18-library example remains pending its separate inputs.

## 5. Are the 18 example libraries available online?

The exact 18 processed libraries are **not confirmed publicly downloadable**. A public Google Drive bundle will be linked here when supplied. Their omission from an annotation or report package alone does not establish whether they exist elsewhere online.

Public ms28 data exist in [GEO GSE279862](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279862), but that series uses a different sample design, including two-replicate anther and male-germ-cell groups. It is not an established match for this tutorial's three-stage, three-replicate collection. Its processed samples were restricted to 18–30 nt, which also does not supply the agreed 18–50-nt test input. See the [GEO sample processing record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8582121).

A public Google Drive bundle can supply the 18 original processed `.tag` files, the matching target sheet, a short processing/provenance README, and SHA-256 checksums. Preserve the existing filenames and file contents. Add its share link here once available; the software repository does not include these libraries. The fixed molecule-oriented Phasis products and candidate FASTAs listed in section 3 remain separate required inputs.
