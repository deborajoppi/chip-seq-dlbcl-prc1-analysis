#!/usr/bin/env bash
set -euo pipefail

mkdir -p results/overlaps

BCOR="results/tmp/normalized_peaks/BCOR_LY1_hg18.bed"
KDM2B="results/tmp/normalized_peaks/KDM2B_LY1_hg18.bed"
H3K27ME3="results/tmp/normalized_peaks/H3K27me3_LY1_hg18.bed"

for file in "${BCOR}" "${KDM2B}" "${H3K27ME3}"; do
  if [[ ! -f "${file}" ]]; then
    printf 'Missing normalized input: %s\nRun make prepare first.\n' "${file}" >&2
    exit 1
  fi
done

bedtools intersect -a "${BCOR}" -b "${KDM2B}" -u > results/overlaps/bcor_kdm2b_overlap.bed
bedtools intersect -a "${BCOR}" -b "${KDM2B}" -v > results/overlaps/bcor_unique_vs_kdm2b.bed
bedtools intersect -a "${KDM2B}" -b "${BCOR}" -v > results/overlaps/kdm2b_unique_vs_bcor.bed
bedtools intersect -a results/overlaps/bcor_kdm2b_overlap.bed -b "${H3K27ME3}" -u > results/overlaps/bcor_kdm2b_overlap_with_h3k27me3.bed
bedtools intersect -a results/overlaps/bcor_kdm2b_overlap.bed -b "${H3K27ME3}" -v > results/overlaps/bcor_kdm2b_overlap_without_h3k27me3.bed

printf 'Peak overlap files written to results/overlaps/\n'
