#!/usr/bin/env bash
set -euo pipefail

mkdir -p results/tmp/normalized_peaks

inputs=(
  "BCOR_LY1_hg18.bed"
  "KDM2B_LY1_hg18.bed"
  "H3K27me3_LY1_hg18.bed"
)

for file in "${inputs[@]}"; do
  src="data/raw/peaks/${file}"
  dst="results/tmp/normalized_peaks/${file}"

  if [[ ! -f "${src}" ]]; then
    printf 'Missing required BED file: %s\n' "${src}" >&2
    exit 1
  fi

  awk 'BEGIN{OFS="\t"} !/^#/ && NF >= 3 {print $1, $2, $3, (NF>=4 ? $4 : "."), (NF>=5 ? $5 : "."), (NF>=6 ? $6 : ".")}' "${src}" \
    | sort -k1,1 -k2,2n -k3,3n > "${dst}"

  printf 'Prepared %s\n' "${dst}"
done
