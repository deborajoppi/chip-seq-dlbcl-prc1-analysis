#!/usr/bin/env bash
set -euo pipefail

# Paths can be overridden from the environment (used by tests/run_tests.sh).
BCOR="${BCOR:-results/tmp/normalized_peaks/BCOR_LY1_hg18.bed}"
KDM2B="${KDM2B:-results/tmp/normalized_peaks/KDM2B_LY1_hg18.bed}"
H3K27ME3="${H3K27ME3:-results/tmp/normalized_peaks/H3K27me3_LY1_hg18.bed}"
OUTDIR="${OUTDIR:-results/overlaps}"

mkdir -p "${OUTDIR}"

for file in "${BCOR}" "${KDM2B}" "${H3K27ME3}"; do
  if [[ ! -f "${file}" ]]; then
    printf 'Missing normalized input: %s\nRun make prepare first.\n' "${file}" >&2
    exit 1
  fi
done

# Peak-level overlap (primary, unit = original peaks):
#   BCOR peaks that touch >=1 KDM2B peak, and vice versa.
bedtools intersect -a "${BCOR}"  -b "${KDM2B}" -u > "${OUTDIR}/bcor_peaks_overlapping_kdm2b.bed"
bedtools intersect -a "${KDM2B}" -b "${BCOR}"  -u > "${OUTDIR}/kdm2b_peaks_overlapping_bcor.bed"
bedtools intersect -a "${BCOR}"  -b "${KDM2B}" -v > "${OUTDIR}/bcor_nonoverlap_strict.bed"
bedtools intersect -a "${KDM2B}" -b "${BCOR}"  -v > "${OUTDIR}/kdm2b_nonoverlap_strict.bed"

# Shared regions (unit = distinct BCOR∩KDM2B intersection fragments).
# This is the definition that reproduces the thesis count of 14,546; it is
# kept because downstream annotation operates on these shared regions, but
# note that fragment counts are not directly comparable to peak counts.
# Inputs are cut to 3 columns first so the column arithmetic below does not
# depend on how many columns the peak files carry.
bedtools intersect -a <(cut -f1-3 "${BCOR}") -b <(cut -f1-3 "${KDM2B}") -wo \
  | awk 'BEGIN{OFS="\t"} {s=($2>$5?$2:$5); e=($3<$6?$3:$6); if (s<e) print $1, s, e, "BCOR_KDM2B_common_region", ".", "."}' \
  | sort -k1,1 -k2,2n -k3,3n \
  | uniq > "${OUTDIR}/bcor_kdm2b_overlap.bed"

bedtools intersect -a "${OUTDIR}/bcor_kdm2b_overlap.bed" -b "${H3K27ME3}" -u > "${OUTDIR}/bcor_kdm2b_overlap_with_h3k27me3.bed"
bedtools intersect -a "${OUTDIR}/bcor_kdm2b_overlap.bed" -b "${H3K27ME3}" -v > "${OUTDIR}/bcor_kdm2b_overlap_without_h3k27me3.bed"

printf 'Peak overlap files written to %s/\n' "${OUTDIR}"
