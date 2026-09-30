#!/usr/bin/env bash
# Copy to your study directory, edit the paths/identities, then source it.
# This file only defines arguments; it does not run an analysis.
SRNA_LIBRARIES=/path/to/processed_libraries
SRNA_TARGETS=/path/to/targets.tsv
SRNA_REFERENCES=/path/to/frozen_reference
SRNA_PHASIS=/path/to/fixed_phasis_products

SRNA_INPUT_ARGS=(
  --processed-libraries "$SRNA_LIBRARIES"
  --library-format tag --min-length 18 --max-length 50
  --targets "$SRNA_TARGETS"
  --assembly-id REPLACE_WITH_ASSEMBLY_ID
  --annotation-id REPLACE_WITH_FROZEN_ANNOTATION_ID
  --reference-orientation molecule
  --reference "miRNA=$SRNA_REFERENCES/miRNA.fa"
  --reference "tRNA=$SRNA_REFERENCES/tRNA.fa"
  --reference "rRNA=$SRNA_REFERENCES/rRNA.fa"
  --reference "snRNA=$SRNA_REFERENCES/snRNA.fa"
  --reference "snoRNA=$SRNA_REFERENCES/snoRNA.fa"
  --phasiRNAs "$SRNA_PHASIS/21_phasiRNAs.tsv"
  --phasiRNAs "$SRNA_PHASIS/24_phasiRNAs.tsv"
)

# Positive log2 fold change means greater abundance in treatment.
# Replace the conditions with names that appear in your targets.tsv.
SRNA_DA_ARGS=(
  --normalization total-qc-depth
  --contrast treatment,control
  --min-count 10 --prefilter-policy all-replicates-in-one-condition --fdr 0.05
)
