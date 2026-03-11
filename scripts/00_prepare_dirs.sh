#!/usr/bin/env bash
set -euo pipefail

mkdir -p \
  data/raw \
  data/raw/peaks \
  data/metadata \
  docs \
  results/annotations \
  results/figures \
  results/overlaps \
  results/tables \
  results/tmp/normalized_peaks \
  reports \
  scripts

printf 'Project directories are ready.\n'
