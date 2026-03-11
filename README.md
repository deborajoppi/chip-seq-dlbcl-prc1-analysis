# BCOR and KDM2B Target Gene Analysis in GCB-DLBCL

This project documents a secondary analysis of public ChIP-seq data from the LY1 cell line to identify putative PRC1.1 target genes in germinal center B-cell-like diffuse large B-cell lymphoma (GCB-DLBCL).

## Overview

Publicly available ChIP-seq peak sets for `KDM2B`, `BCOR`, and `H3K27me3` were used to define chromatin regions bound by PRC1.1-associated factors while excluding regions likely shared with PRC2.

The analysis logic is:

1. start from public LY1 ChIP-seq peak data aligned to `hg18`
2. annotate peaks with `HOMER`
3. identify overlapping `KDM2B` and `BCOR` peaks
4. use `H3K27me3` to flag PRC2-associated regions for exclusion
5. summarize genomic annotations and enriched biological processes for the retained gene set

For the BCOR-KDM2B overlap, the rerun workflow uses distinct genomic intersection regions between the two public peak sets. With the currently downloaded public files, that reproduces the thesis-reported overlap count of `14,546`.

## Public Data Used

The public datasets used in the thesis analysis are:

| Assay | Target | GEO accession | URL | Notes |
|---|---|---|---|---|
| ChIP-seq | BCOR | `GSE29282` | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29282 | Source used for BCOR LY1 peaks |
| ChIP-seq | KDM2B | `GSE81623` | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81623 | Source used for KDM2B LY1 peaks |
| ChIP-seq | H3K27me3 | `GSM763414` | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM763414 | Source used to exclude PRC2-associated regions |

Because GEO releases can package peak files differently across accessions, this repo standardizes the rerun workflow around three local BED inputs:

- `data/raw/peaks/BCOR_LY1_hg18.bed`
- `data/raw/peaks/KDM2B_LY1_hg18.bed`
- `data/raw/peaks/H3K27me3_LY1_hg18.bed`

The step-by-step workflow for preparing those files and rerunning the analysis is documented in `docs/reanalysis_workflow.md`.

## Key Results

| Target set | Total peaks |
|---|---:|
| BCOR | 17,549 |
| KDM2B | 26,471 |
| BCOR and KDM2B overlap | 14,546 |
| H3K27me3 | 39,158,099 |

Additional reported findings:

- `14,546` peaks overlapped between `BCOR` and `KDM2B`.
- `6,943` overlapping peaks were located in promoter-TSS and intronic regions.
- `11,925` peaks were unique to `KDM2B`.
- `3,002` peaks were unique to `BCOR`.
- Enriched biological processes included proteasome-mediated ubiquitin-dependent protein catabolic process, mitotic cell cycle phase transition, Wnt signaling pathway, nucleocytoplasmic transport, and intrinsic apoptotic signaling pathway.

Representative genes highlighted in the thesis text include `PRICKLE1`, `CUL1`, `MARCHF6`, `FBXL3`, `SMAD7`, `CLOCK`, `POU4F1`, `STK24`, `MYC`, and `JAK2`.

Current public rerun status:

- the overlap count `14,546` is reproduced exactly from the public files
- the Venn-style unique counts are also reproduced (`3,002` BCOR-unique and `11,925` KDM2B-unique)
- the promoter-plus-intron figure from the thesis text is not yet reproduced by the current public rerun annotation workflow
- `H3K27me3` threshold checks from `0.1` to `0.5` only modestly change the retained overlap-without-H3K27me3 count, so the unresolved `6,943` value likely reflects a different filtering or reporting choice

## Repository Layout

```text
chip-seq-dlbcl-prc1-analysis/
├── docs/
│   └── reanalysis_workflow.md
├── data/
│   ├── metadata/
│   │   └── public_sources.csv
│   └── raw/
│       ├── .gitkeep
│       └── peaks/
│           └── .gitkeep
├── reports/
│   └── report.Rmd
├── results/
│   ├── annotations/
│   │   └── .gitkeep
│   ├── figures/
│   │   └── .gitkeep
│   ├── overlaps/
│   │   └── .gitkeep
│   ├── tables/
│   │   ├── candidate_target_genes.csv
│   │   ├── genomic_region_summary.csv
│   │   ├── go_terms_summary.csv
│   │   ├── peak_totals.csv
│   │   └── thesis_count_comparison.csv
│   └── tmp/
│       └── .gitkeep
├── scripts/
│   ├── 00_prepare_dirs.sh
│   ├── 00_fetch_public_data.sh
│   ├── 01_prepare_peak_inputs.sh
│   ├── 02_intersect_peaks.sh
│   ├── 03_annotate_peaks.sh
│   ├── 04_summarize_annotations.py
│   ├── 05_run_go_enrichment.R
│   ├── 06_render_report.R
│   ├── 07_h3k27me3_threshold_sweep.py
│   ├── 08_make_figures.R
│   └── wig_to_bed.py
├── .gitignore
├── LICENSE
├── Makefile
└── README.md
```

## Run

Prepare directories:

```bash
make init
```

Run the BED-based reanalysis after placing the three peak files in `data/raw/peaks/`:

```bash
make fetch
make prepare
make overlap
make annotate
make summarize
make figures
make sweep
make report
```

Optional GO enrichment:

```bash
make go
```

The rendered output is written to `results/report.html`. Detailed rerun notes are in `docs/reanalysis_workflow.md`.

## Data Notes

- This repository is intended for `public` data and derived summaries only.
- Raw downloaded files are excluded from git by default under `data/raw/`.
- GEO accessions and source links are tracked in `data/metadata/public_sources.csv`.
- The current committed summary tables capture the thesis-reported values and provide expected targets for rerunning the pipeline.
- The workflow assumes that the public peak files used for rerun are aligned to `hg18`, matching the thesis analysis.
- `H3K27me3` is distributed publicly as a WIG signal track, so the exclusion BED used in the rerun is derived from signal intervals rather than a deposited peak BED.
- The overlap count is computed from distinct BCOR-KDM2B intersection regions. Strict non-overlap peak counts are tracked separately because they are not identical to the Venn-style counts reported in the thesis.
- `make figures` generates PNG plots in `results/figures/` for the rerun counts, region distribution, and H3K27me3 threshold sensitivity.
- `make figures` also generates a top-candidate-gene bar plot and a gene-by-set heatmap from the shared non-H3K27me3 annotation output.
- It also generates a promoter-restricted top-target plot excluding the self-genes `BCOR` and `KDM2B`.
- `make gsea` runs a preranked Hallmark GSEA on promoter-focused shared non-H3K27me3 target genes and writes a pathway plot plus result table.
- `make sweep` runs a small sensitivity analysis over multiple `H3K27me3` signal thresholds and writes `results/tables/h3k27me3_threshold_sweep.csv`.
