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
7. testing CpG-island enrichment and running GO/Hallmark over-representation against explicit backgrounds

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
- R packages `clusterProfiler`, `org.Hs.eg.db`, `AnnotationDbi`, `msigdbr`, `rmarkdown`

All of these are pinned in `environment.yml`.

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

- `results/overlaps/bcor_peaks_overlapping_kdm2b.bed` (peak level)
- `results/overlaps/kdm2b_peaks_overlapping_bcor.bed` (peak level)
- `results/overlaps/bcor_kdm2b_overlap.bed` (shared intersection regions)
- `results/overlaps/bcor_nonoverlap_strict.bed`
- `results/overlaps/kdm2b_nonoverlap_strict.bed`
- `results/overlaps/bcor_kdm2b_overlap_with_h3k27me3.bed`
- `results/overlaps/bcor_kdm2b_overlap_without_h3k27me3.bed`

Two overlap definitions are produced. Peak-level files count original peaks that touch at least one peak of the other factor (14,118 BCOR, 12,667 KDM2B). The shared-region file contains distinct BCOR∩KDM2B intersection fragments (14,546, matching the thesis) and is what gets annotated downstream. The two are in different units and should not be subtracted from each other.

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

### 6. CpG-island background and enrichment

```bash
make cgi
make enrich
```

`make cgi` measures the fraction of each peak set on CpG islands against 100 genome shuffles (same chromosome, same interval lengths) and builds a gene-level promoter table (`results/tables/promoter_universe.tsv`) marking CGI promoters and bound promoters.

`make enrich` runs GO BP and Hallmark over-representation for promoters bound by shared regions outside the H3K27me3 mask, against all genes and against CGI-promoter genes only. hg18 gene symbols are mapped to current Entrez IDs through RefSeq accessions.

### 7. Sweep H3K27me3 thresholds

```bash
make sweep
```

This writes `results/tables/h3k27me3_threshold_sweep.csv` and helps quantify how sensitive the PRC2-exclusion step is to the WIG-to-BED threshold.

Current interpretation from the public rerun:

- the shared BCOR-KDM2B overlap count is reproduced exactly at `14,546`
- the Venn-style unique counts are also reproduced (`3,002` BCOR-unique and `11,925` KDM2B-unique)
- the promoter-plus-intron figure from the thesis text is **not** currently reproduced by the public rerun annotation workflow
- threshold sweeps from `0.1` to `0.5` change the retained count only modestly, but at `1.0` it nearly doubles (6,714 → 11,475), so the H3K27me3 filter is threshold-dependent; the `6,943` thesis value is not explained by any threshold tested

### 8. Render the report

```bash
make report
```

This writes `results/report.html`.

## Interpretation Notes

- The thesis summary reports `14,546` overlapping peaks between `BCOR` and `KDM2B`.
- It also reports `11,925` `KDM2B`-unique and `3,002` `BCOR`-unique peaks, computed as total peaks minus shared regions. Strict per-peak counts are `13,804` and `3,430`.
- `H3K27me3` is used here as an exclusion filter for candidate regions likely co-occupied by `PRC2`.

If your rerun counts differ, record the reason in the README:

- a different supplementary file was used
- the public files were not the same processed peak set as in the thesis
- chromosome naming required normalization
- overlap logic differed from the original run
- the original thesis analysis excluded additional categories not yet scripted here
