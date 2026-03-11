# BCOR and KDM2B Target Gene Analysis in GCB-DLBCL

Reproducible repository scaffold for the thesis analysis section:

`1.19 Identification of BCOR and KDM2B target genes in GCB-type DLBCL`

This project documents a secondary analysis of public ChIP-seq data from the LY1 cell line to identify putative PRC1.1 target genes in germinal center B-cell-like diffuse large B-cell lymphoma (GCB-DLBCL).

## Overview

Publicly available ChIP-seq peak sets for `KDM2B`, `BCOR`, and `H3K27me3` were used to define chromatin regions bound by PRC1.1-associated factors while excluding regions likely shared with PRC2.

The analysis logic is:

1. start from public LY1 ChIP-seq peak data aligned to `hg18`
2. annotate peaks with `HOMER`
3. identify overlapping `KDM2B` and `BCOR` peaks
4. use `H3K27me3` to flag PRC2-associated regions for exclusion
5. summarize genomic annotations and enriched biological processes for the retained gene set

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

## Repository Layout

```text
chip-seq-dlbcl-prc1-analysis/
├── data/
│   ├── metadata/
│   │   └── public_sources_template.csv
│   └── raw/
│       └── .gitkeep
├── reports/
│   └── report.Rmd
├── results/
│   ├── figures/
│   │   └── .gitkeep
│   ├── tables/
│   │   ├── go_terms_summary.csv
│   │   └── peak_totals.csv
│   └── tmp/
│       └── .gitkeep
├── scripts/
│   ├── 00_prepare_dirs.sh
│   └── 01_render_report.R
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

Render the summary report:

```bash
make report
```

The rendered output is written to `results/report.html`.

## Data Notes

- This repository is intended for `public` data and derived summaries only.
- Raw downloaded files are excluded from git by default under `data/raw/`.
- The current committed tables capture the summary values reported in the thesis text. You can replace or extend them once the original analysis inputs and figure panels are added.

## Next Additions

- add the original public ChIP-seq source accessions and download links to `data/metadata/public_sources_template.csv`
- add the figure panel files corresponding to thesis Figure 31 into `results/figures/`
- add the original overlap tables, genomic annotation tables, and GO output tables if available
- expand the report to include the exact filtering rules used to exclude `H3K27me3`-associated peaks
