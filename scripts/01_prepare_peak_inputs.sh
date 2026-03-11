#!/usr/bin/env bash
set -euo pipefail

mkdir -p results/tmp/normalized_peaks
mkdir -p data/raw/peaks

BCOR_BED="data/raw/peaks/BCOR_LY1_hg18.bed"
KDM2B_BED="data/raw/peaks/KDM2B_LY1_hg18.bed"
H3K27ME3_BED="data/raw/peaks/H3K27me3_LY1_hg18.bed"

if [[ ! -f "${BCOR_BED}" ]]; then
  if [[ ! -f data/raw/peaks/BCOR_LY1_hg18.source.bed.gz ]]; then
    printf 'Missing BCOR source file. Run make fetch or place %s manually.\n' "${BCOR_BED}" >&2
    exit 1
  fi
  gunzip -c data/raw/peaks/BCOR_LY1_hg18.source.bed.gz \
    | awk 'BEGIN{OFS="\t"} $1 !~ /^track/ && NF >= 3 {print $1, $2, $3, "BCOR_peak_" NR, ".", "."}' \
    > "${BCOR_BED}"
fi

if [[ ! -f "${KDM2B_BED}" ]]; then
  if [[ ! -f data/raw/peaks/KDM2B_LY1_hg18.source.txt.gz ]]; then
    printf 'Missing KDM2B source file. Run make fetch or place %s manually.\n' "${KDM2B_BED}" >&2
    exit 1
  fi
  gunzip -c data/raw/peaks/KDM2B_LY1_hg18.source.txt.gz \
    | awk 'BEGIN{OFS="\t"} NF >= 3 {print $1, $2, $3, "KDM2B_peak_" NR, $5, "."}' \
    > "${KDM2B_BED}"
fi

if [[ ! -f "${H3K27ME3_BED}" ]]; then
  if [[ ! -f data/raw/peaks/H3K27me3_LY1_hg18.source.wig.gz ]]; then
    printf 'Missing H3K27me3 source file. Run make fetch or place %s manually.\n' "${H3K27ME3_BED}" >&2
    exit 1
  fi
  python3 scripts/wig_to_bed.py \
    --input data/raw/peaks/H3K27me3_LY1_hg18.source.wig.gz \
    --output "${H3K27ME3_BED}" \
    --min-signal "${MIN_H3K27ME3_SIGNAL:-0.1}"
fi

inputs=("${BCOR_BED}" "${KDM2B_BED}" "${H3K27ME3_BED}")

for src in "${inputs[@]}"; do
  file="$(basename "${src}")"
  dst="results/tmp/normalized_peaks/${file}"
  awk 'BEGIN{OFS="\t"} !/^#/ && NF >= 3 {print $1, $2, $3, (NF>=4 ? $4 : "."), (NF>=5 ? $5 : "."), (NF>=6 ? $6 : ".")}' "${src}" \
    | sort -k1,1 -k2,2n -k3,3n > "${dst}"
  printf 'Prepared %s\n' "${dst}"
done
