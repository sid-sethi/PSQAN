# PSQAN - Post Sqanti QC ANalysis

<!-- badges: start -->
![Maintainer](https://img.shields.io/badge/maintainer-SidSethi-blue)
![Generic badge](https://img.shields.io/badge/WMS-snakemake-blue.svg)
![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)
![Linux](https://svgshare.com/i/Zhy.svg)
![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
[![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/sid-sethi/APTARS/blob/main/LICENSE)
<!-- badges: end -->

PSQAN is a `snakemake` pipeline for performing post Sqanti QC analysis on Pacbio sequencing data.

<p align="center">
  <img src="dag/dag.png" width="500" height="300"/>  
</p>


# Getting Started

## Input

- Sqanti _classification.txt file
- Sqanti _corrected.fasta file
- Gene of interest (optional)


## Depedencies

- [miniconda](https://conda.io/miniconda.html)
- The rest of the dependencies (including snakemake) are installed via conda through the `environment.yml` file

## Installation

Clone the directory:

```bash
git clone --recursive https://github.com/sid-sethi/PSQAN.git
```

Create conda environment for the pipeline which will install all the dependencies:

```bash
cd PSQAN
conda env create -f environment.yml
```

## Usage

Edit `config.yml` to set up the working directory and input files/directories. `snakemake` command should be issued from within the pipeline directory. Please note that before you run any of the `snakemake` commands, make sure to first activate the conda environment using the command `conda activate psqan`.

```bash
cd PSQAN
conda activate psqan
snakemake --use-conda -j <num_cores> all
```
It is a good idea to do a dry run (using -n parameter) to view what would be done by the pipeline before executing the pipeline.

```bash
snakemake --use-conda -n all
```

You can visualise the processes to be executed in a DAG:

```bash
snakemake --dag | dot -Tpng > dag.png
```

To exit a running `snakemake` pipeline, hit `ctrl+c` on the terminal. If the pipeline is running in the background, you can send a `TERM` signal which will stop the scheduling of new jobs and wait for all running jobs to be finished.

```bash
killall -TERM snakemake
```

To deactivate the conda environment:
```bash
conda deactivate
```
