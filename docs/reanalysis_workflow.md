# Reanalysis Workflow

This document turns the thesis subsection into a rerunnable workflow that can be shown step by step on GitHub.

## Goal

Identify putative `PRC1.1` target genes in `LY1` cells by:

1. taking public `BCOR` and `KDM2B` ChIP-seq peak sets
2. finding the shared peaks between those factors
3. identifying which shared peaks also overlap `H3K27me3`
4. retaining the shared non-`H3K27me3` peaks as PRC1.1-focused candidates
5. annotating peaks to nearby genes with `HOMER`
6. summarizing genomic distribution and candidate target genes
7. optionally running GO enrichment on the candidate gene list

## Public Sources

The accessions used in the thesis analysis are recorded in `data/metadata/public_sources.csv`:

- `GSE29282` for `BCOR`
- `GSE81623` for `KDM2B`
- `GSM763414` for `H3K27me3`

## Input Convention

The repo can download the public GEO source files directly:

```bash
make fetch
```

This retrieves:

- the `BCOR` overlap peak BED-like file
- the `KDM2B` peak table
- the `H3K27me3` read-density WIG

If you prefer to provide files manually, place BED files in `data/raw/peaks/` using these filenames:

- `BCOR_LY1_hg18.bed`
- `KDM2B_LY1_hg18.bed`
- `H3K27me3_LY1_hg18.bed`

Expected BED columns:

1. chromosome
2. start
3. end
4. optional peak name
5. optional score
6. optional strand

Only the first three columns are required by the overlap steps.

## Tool Requirements

- `bedtools`
- `annotatePeaks.pl` from HOMER
- `python3`
- `Rscript`
- optional for GO: `clusterProfiler` and `org.Hs.eg.db`

## Step-by-Step Commands

### 1. Prepare the repository structure

```bash
make init
```

### 2. Normalize and sort the BED inputs

```bash
make prepare
```

This step:

- converts the downloaded public source files into standard BED inputs if needed
- keeps the first six BED columns if extra columns are present
- sorts each file by chromosome and genomic coordinate
- writes standardized copies to `results/tmp/normalized_peaks/`

For `H3K27me3`, GEO exposes a `variableStep` WIG signal track rather than a called-peak BED. The current workflow converts that signal to merged exclusion intervals using a configurable minimum signal threshold:

```bash
MIN_H3K27ME3_SIGNAL=0.1 make prepare
```

This is an approximation of the exclusion mask and may not exactly match the original thesis run if the original analysis used a different peak-calling or thresholding strategy.

### 3. Identify overlaps

```bash
make overlap
```

This step creates:

- `results/overlaps/bcor_kdm2b_overlap.bed`
- `results/overlaps/bcor_nonoverlap_strict.bed`
- `results/overlaps/kdm2b_nonoverlap_strict.bed`
- `results/overlaps/bcor_kdm2b_overlap_with_h3k27me3.bed`
- `results/overlaps/bcor_kdm2b_overlap_without_h3k27me3.bed`

The primary overlap file is built from distinct BCOR-KDM2B genomic intersection regions, not just from an asymmetric `bedtools intersect -u` call. With the public files used here, that reproduces the thesis overlap count of `14,546`.

### 4. Annotate peak sets with HOMER

```bash
make annotate
```

This step annotates:

- `BCOR`
- `KDM2B`
- `BCOR-KDM2B overlap`
- `BCOR-KDM2B overlap without H3K27me3`

Outputs are written to `results/annotations/`.

### 5. Build summary tables

```bash
make summarize
```

This creates:

- `results/tables/peak_totals.csv`
- `results/tables/genomic_region_summary.csv`
- `results/tables/candidate_target_genes.csv`
- `results/tables/thesis_count_comparison.csv`

### 6. Optional: GO enrichment

```bash
make go
```

If the required R packages are installed, this reads `results/tables/candidate_target_genes.csv` and writes GO enrichment results to `results/tables/go_enrichment.csv`.

### 7. Render the report

```bash
make report
```

This writes `results/report.html`.

## Interpretation Notes

- The thesis summary reports `14,546` overlapping peaks between `BCOR` and `KDM2B`.
- It also reports `11,925` `KDM2B`-unique peaks and `3,002` `BCOR`-unique peaks.
- `H3K27me3` is used here as an exclusion filter for candidate regions likely co-occupied by `PRC2`.

If your rerun counts differ, record the reason in the README:

- a different supplementary file was used
- the public files were not the same processed peak set as in the thesis
- chromosome naming required normalization
- overlap logic differed from the original run
- the original thesis analysis excluded additional categories not yet scripted here
