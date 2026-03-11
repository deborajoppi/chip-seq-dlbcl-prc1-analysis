#!/usr/bin/env bash
set -euo pipefail

mkdir -p results/annotations

GENOME="${1:-hg18}"

peak_sets=(
  "BCOR:results/tmp/normalized_peaks/BCOR_LY1_hg18.bed"
  "KDM2B:results/tmp/normalized_peaks/KDM2B_LY1_hg18.bed"
  "BCOR_KDM2B_overlap:results/overlaps/bcor_kdm2b_overlap.bed"
  "BCOR_KDM2B_overlap_no_H3K27me3:results/overlaps/bcor_kdm2b_overlap_without_h3k27me3.bed"
)

for entry in "${peak_sets[@]}"; do
  label="${entry%%:*}"
  input="${entry#*:}"
  output="results/annotations/${label}.annotated.tsv"

  if [[ ! -f "${input}" ]]; then
    printf 'Missing peak file for annotation: %s\n' "${input}" >&2
    exit 1
  fi

  annotatePeaks.pl "${input}" "${GENOME}" > "${output}"
  printf 'Annotated %s\n' "${output}"
done
