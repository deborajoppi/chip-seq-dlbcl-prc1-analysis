SHELL := /bin/bash

.PHONY: init report

init:
	bash scripts/00_prepare_dirs.sh

report:
	Rscript scripts/01_render_report.R
