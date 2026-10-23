SHELL := /bin/bash
.DELETE_ON_ERROR:

RAW   := data/raw/peaks
ANN   := data/raw/annotation
NORM  := results/tmp/normalized_peaks
OVL   := results/overlaps
TAB   := results/tables

SOURCES    := $(RAW)/BCOR_LY1_hg18.source.bed.gz $(RAW)/KDM2B_LY1_hg18.source.txt.gz $(RAW)/H3K27me3_LY1_hg18.source.wig.gz
NORMALIZED := $(NORM)/BCOR_LY1_hg18.bed $(NORM)/KDM2B_LY1_hg18.bed $(NORM)/H3K27me3_LY1_hg18.bed
OVERLAPS   := $(OVL)/bcor_kdm2b_overlap.bed $(OVL)/bcor_kdm2b_overlap_without_h3k27me3.bed
ANNOTATED  := results/annotations/BCOR_KDM2B_overlap_no_H3K27me3.annotated.tsv
SUMMARY    := $(TAB)/peak_totals.csv
SWEEP      := $(TAB)/h3k27me3_threshold_sweep.csv
CGI        := $(TAB)/promoter_universe.tsv
ENRICH     := $(TAB)/enrichment_summary.csv
REPORT     := results/report.html

.PHONY: all init fetch prepare overlap annotate summarize sweep cgi enrich figures report test clean

all: $(REPORT)

init:
	bash scripts/00_prepare_dirs.sh

fetch: $(SOURCES)
$(SOURCES) &:
	bash scripts/00_fetch_public_data.sh

prepare: $(NORMALIZED)
$(NORMALIZED) &: $(SOURCES) scripts/01_prepare_peak_inputs.sh scripts/wig_to_bed.py
	bash scripts/01_prepare_peak_inputs.sh

overlap: $(OVERLAPS)
$(OVERLAPS) &: $(NORMALIZED) scripts/02_intersect_peaks.sh
	bash scripts/02_intersect_peaks.sh

annotate: $(ANNOTATED)
$(ANNOTATED): $(OVERLAPS) scripts/03_annotate_peaks.sh
	bash scripts/03_annotate_peaks.sh

summarize: $(SUMMARY)
$(SUMMARY): $(ANNOTATED) scripts/04_summarize_annotations.py
	python3 scripts/04_summarize_annotations.py

sweep: $(SWEEP)
$(SWEEP): $(OVERLAPS) scripts/07_h3k27me3_threshold_sweep.py
	python3 scripts/07_h3k27me3_threshold_sweep.py

cgi: $(CGI)
$(CGI): $(OVERLAPS) scripts/10_cgi_background.py
	python3 scripts/10_cgi_background.py

enrich: $(ENRICH)
$(ENRICH): $(CGI) scripts/05_run_enrichment.R
	Rscript scripts/05_run_enrichment.R

figures: $(SUMMARY) $(SWEEP)
	Rscript scripts/08_make_figures.R

$(REPORT): $(SUMMARY) $(SWEEP) $(CGI) $(ENRICH) reports/report.Rmd
	Rscript scripts/08_make_figures.R
	Rscript scripts/06_render_report.R

report: $(REPORT)

test:
	bash tests/run_tests.sh

clean:
	rm -rf results/tmp/* $(OVL)/*.bed
