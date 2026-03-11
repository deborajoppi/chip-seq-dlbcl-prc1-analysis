#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/raw/peaks

curl -sSL -o data/raw/peaks/BCOR_LY1_hg18.source.bed.gz \
  'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM763nnn/GSM763410/suppl/GSM763410_OCI-LY1_BCOR_peaks_overlap.R1-R2.txt.wgl.bed.gz'

curl -sSL -o data/raw/peaks/KDM2B_LY1_hg18.source.txt.gz \
  'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2171nnn/GSM2171650/suppl/GSM2171650_Sample_1_KDM2B_ChIPseeqer_T15F2peaks.txt.gz'

curl -sSL -o data/raw/peaks/H3K27me3_LY1_hg18.source.wig.gz \
  'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM763nnn/GSM763414/suppl/GSM763414_OCI-LY1_H3K27me3_CHIP_ReadDens.wig.gz'

printf 'Downloaded public GEO source files into data/raw/peaks/.\n'
