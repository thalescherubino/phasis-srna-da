# Download the B73v5 reference bundle

The tutorial needs one small, fixed reference bundle. It contains the B73v5
small-RNA annotation FASTA files and the mandatory Phasis phasiRNA tables. It
does not contain sequencing libraries.

Use this guide after [setting up the software](phasis_srna_da_setup.md), in an
existing `srun` compute shell.

## Download and check the bundle

When the bundle has been uploaded, replace the two values below with its public
release URL and SHA-256 checksum.

```bash
export SRNA_HOME="/quobyte/bcmeyersgrp/$USER/srna-da"
export REFERENCE_URL="PASTE_THE_RELEASE_ASSET_URL_HERE"
export REFERENCE_SHA256="PASTE_THE_ARCHIVE_SHA256_HERE"

mkdir -p "$SRNA_HOME"
curl -L "$REFERENCE_URL" -o "$SRNA_HOME/ms28-example-reference-v1.zip"
echo "$REFERENCE_SHA256  $SRNA_HOME/ms28-example-reference-v1.zip" | sha256sum --check
unzip "$SRNA_HOME/ms28-example-reference-v1.zip" -d "$SRNA_HOME"
sha256sum --check "$SRNA_HOME/ms28-example/SHA256SUMS.references"
```

After extraction, the directory must look like this:

```text
ms28-example/
├── annotation/fasta/
├── candidate_fasta/
├── phasiRNAs/
├── README.references.md
└── SHA256SUMS.references
```

The final CLI configuration uses `SRNA_EXAMPLE_ROOT` set to
`$SRNA_HOME/ms28-example`. A separate library archive can be extracted to the
same directory later; it adds `libraries/` and `targets.tsv` without replacing
the reference files.

## Contents required in the reference bundle

Package exactly these local inputs, retaining the filenames shown:

```text
annotation/fasta/miRNA_bona_fide.mature_and_star.fa
annotation/fasta/tRNA.CCA_aware_quantification.fa
annotation/fasta/rRNA.fa
annotation/fasta/snRNA.fa
annotation/fasta/snoRNA.fa
annotation/fasta/miRNA_unsupported_candidates.fa
candidate_fasta/21_PHAS_like.molecule_oriented.fa
candidate_fasta/24_PHAS_like.molecule_oriented.fa
phasiRNAs/21_PHAS.molecule_oriented.phasiRNAs.tsv
phasiRNAs/24_PHAS.molecule_oriented.phasiRNAs.tsv
```

Include `README.references.md` and `SHA256SUMS.references` in the archive. The
README must identify the frozen assembly as `Zm-B73-REFERENCE-NAM-5.0` and the
reference set as
`B73v5_curated_sRNA_reference_v1_2026-08-27__PHAS_molecule_oriented_v1_2026-08-30`.

The ordinary PHAS FASTA files, GFF3 files, plots, discovery outputs, and other
reference-package material are not inputs to this differential-abundance
tutorial and should stay out of this bundle.

## Where to put the files

Publish the reference ZIP as a GitHub Release asset in the public
`thalescherubino/phasis-srna-da` repository. It contains reference inputs only
and should be accompanied by its archive SHA-256 checksum.

The sequencing-library archive is separate and must remain private. A private
Google Drive folder is the simplest option. If encryption is required, create
one AES-256 encrypted `.7z` archive containing:

```text
ms28-example/
├── libraries/            # the 18 processed tag libraries
├── targets.tsv
├── README.libraries.md
└── SHA256SUMS.libraries
```

Share the password through a password manager or another secure channel; never
place it in GitHub, the archive filename, or this tutorial. The next analyst
needs the private download link, any access instructions, the password, and the
archive checksum.
