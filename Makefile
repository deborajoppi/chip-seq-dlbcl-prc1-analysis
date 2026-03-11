SHELL := /bin/bash

.PHONY: init fetch prepare overlap annotate summarize go report all

init:
	bash scripts/00_prepare_dirs.sh

fetch:
	bash scripts/00_fetch_public_data.sh

prepare:
	bash scripts/01_prepare_peak_inputs.sh

overlap:
	bash scripts/02_intersect_peaks.sh

annotate:
	bash scripts/03_annotate_peaks.sh

summarize:
	python3 scripts/04_summarize_annotations.py

go:
	Rscript scripts/05_run_go_enrichment.R

report:
	Rscript scripts/06_render_report.R

all: init fetch prepare overlap annotate summarize report
