#!/usr/bin/env bash
# Copy to your study directory and set the lab-provided example location.
# Input files are read-only; this file does not execute any analysis.
SRNA_EXAMPLE_ROOT=/path/to/lab/ms28-example
SRNA_REFERENCE="$SRNA_EXAMPLE_ROOT/annotation/fasta"
SRNA_PHAS="$SRNA_EXAMPLE_ROOT/phasiRNAs"
SRNA_CANDIDATES="$SRNA_EXAMPLE_ROOT/candidate_fasta"

SRNA_INPUT_ARGS=(
  --processed-libraries "$SRNA_EXAMPLE_ROOT/libraries"
  --library-format tag --min-length 18 --max-length 50
  --targets "$SRNA_EXAMPLE_ROOT/targets.tsv"
  --assembly-id Zm-B73-REFERENCE-NAM-5.0
  --annotation-id B73v5_curated_sRNA_reference_v1_2026-08-27__PHAS_molecule_oriented_v1_2026-08-30
  --reference-orientation molecule
  --reference "miRNA_bona_fide=$SRNA_REFERENCE/miRNA_bona_fide.mature_and_star.fa"
  --reference "tRNA=$SRNA_REFERENCE/tRNA.CCA_aware_quantification.fa"
  --reference "rRNA=$SRNA_REFERENCE/rRNA.fa"
  --reference "snRNA=$SRNA_REFERENCE/snRNA.fa"
  --reference "snoRNA=$SRNA_REFERENCE/snoRNA.fa"
  --reference "unsupported_MIR_candidate=$SRNA_REFERENCE/miRNA_unsupported_candidates.fa"
  --reference "PHAS_like-21=$SRNA_CANDIDATES/21_PHAS_like.molecule_oriented.fa"
  --reference "PHAS_like-24=$SRNA_CANDIDATES/24_PHAS_like.molecule_oriented.fa"
  --phasiRNAs "$SRNA_PHAS/21_PHAS.molecule_oriented.phasiRNAs.tsv"
  --phasiRNAs "$SRNA_PHAS/24_PHAS.molecule_oriented.phasiRNAs.tsv"
)

SRNA_DA_ARGS=(
  --normalization total-qc-depth
  --contrast ms28_0p4,WT_0p4 --contrast ms28_2,WT_2 --contrast ms28_5,WT_5
  --contrast WT_2,WT_0p4 --contrast WT_5,WT_0p4 --contrast WT_5,WT_2
  --contrast ms28_2,ms28_0p4 --contrast ms28_5,ms28_0p4 --contrast ms28_5,ms28_2
  --min-count 10 --prefilter-policy all-replicates-in-one-condition --fdr 0.05
)
