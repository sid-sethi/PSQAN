## -------------------------------------------------------------------------
## Makefile to help run PSQAN
## -------------------------------------------------------------------------

SHELL=/bin/bash
CONDA_BASE=$$(conda info --base)
CONDA_ACTIVATE=source ${CONDA_BASE}/etc/profile.d/conda.sh ; conda activate ; conda activate

help:  ## Show this help message
	@sed -ne '/@sed/!s/## //p' $(MAKEFILE_LIST)

run_psqan: ## run PSQAN
	@($(CONDA_ACTIVATE) ./psqan_venv;\
	snakemake -j 10 all)

run_report: ## generate snakemake reort
	@($(CONDA_ACTIVATE) ./psqan_venv;\
	snakemake --report report.html)

unit_tests: ## run all unit tests
	@($(CONDA_ACTIVATE) ./psqan_venv;\
	Rscript tests/testthat.R)