#!/usr/bin/env bash
set -euo pipefail

# ChIP-seq inputs from GEO
mkdir -p data/raw/peaks data/raw/annotation

curl -sSL -o data/raw/peaks/BCOR_LY1_hg18.source.bed.gz \
  'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM763nnn/GSM763410/suppl/GSM763410_OCI-LY1_BCOR_peaks_overlap.R1-R2.txt.wgl.bed.gz'

curl -sSL -o data/raw/peaks/KDM2B_LY1_hg18.source.txt.gz \
  'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2171nnn/GSM2171650/suppl/GSM2171650_Sample_1_KDM2B_ChIPseeqer_T15F2peaks.txt.gz'

curl -sSL -o data/raw/peaks/H3K27me3_LY1_hg18.source.wig.gz \
  'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM763nnn/GSM763414/suppl/GSM763414_OCI-LY1_H3K27me3_CHIP_ReadDens.wig.gz'

# hg18 reference annotation from UCSC, used for the CpG-island background
# analysis (scripts/10_cgi_background.py) and for genome-fraction reporting
# in the H3K27me3 threshold sweep.
UCSC='https://hgdownload.soe.ucsc.edu/goldenPath/hg18/database'
curl -sSL -o data/raw/annotation/hg18_chromInfo.txt.gz   "${UCSC}/chromInfo.txt.gz"
curl -sSL -o data/raw/annotation/hg18_cpgIslandExt.txt.gz "${UCSC}/cpgIslandExt.txt.gz"
curl -sSL -o data/raw/annotation/hg18_refGene.txt.gz     "${UCSC}/refGene.txt.gz"

printf 'Downloaded GEO peak files into data/raw/peaks/ and UCSC hg18 annotation into data/raw/annotation/.\n'
