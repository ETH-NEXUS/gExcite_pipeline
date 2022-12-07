[![Snakemake](https://img.shields.io/badge/snakemake-≥6.12.1-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![Code style: snakefmt](https://img.shields.io/badge/code%20style-snakefmt-000000.svg)](https://github.com/snakemake/snakefmt)
[![GitHub Super-Linter](https://github.com/ETH-NEXUS/scGATE_workflow/workflows/Lint%20Code%20Base/badge.svg)](https://github.com/marketplace/actions/super-linter)
# scGate 
![Workflow Figure](https://github.com/ETH-NEXUS/gExcite_pipeline/blob/main/images/Workflow_Figure_gExcite.png)

## General overview
gExcite is a start-to-end workflow embedded in Snakemake that provides both, gene expression and CITE-seq analysis, as well as hashing deconvolution.  
For an overview of all steps please see the Snakemake [rulegraph](images/gExcite_pipeline_rulegraph.png).

## Installation instructions
### Pipeline
Given conda is installed on your system the pipeline can be set up using `snakedeploy`.

First, create and activate an environment including Mamba, Snakemake and Snakedeploy:

```
conda create -c bioconda -c conda-forge --name snakemake mamba snakemake snakedeploy
conda activate snakemake
```

Snakedeploy can now be used to deploy the workflow: 

```
snakedeploy deploy-workflow https://github.com/ETH-NEXUSgExcite_pipeline --tag main .
```

### Dependencies
Most of the software used in the default workflow can be installed in an automated fashion using Snakemake's [--use-conda](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management) functionality. 
The following software needs to be installed manually.

- [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger): Follow the instructions on the 10xGenomics installation support page to install cellranger and to include the cellranger binary to your path.
Webpage: [https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation)


## Example data
We provide example data [here](insert link) that can be used for a test run.  
It consists of three hashed samples of human PBMC cells.


## Before running the pipeline
Before the pipeline can be run make sure that
- ADT and the GEX FASTQ files are provided in the folder structure specified below
- the pipeline is configured to your data
- necessary preprocessing was performed

### Prepare the FASTQ files
The pipeline expects the FASTQ files per sample to be in the following folder structure and adhere to the naming schema: 

/path/to/**input_fastqs_gex**/**SAMPLENAME**/**SAMPLENAME**\_S[Number]\_L00[Lane Number]\_[Read Type]\_001.fastq.gz
```
input_fastqs_adt
└── SAMPLENAME
    ├── SAMPLENAME_S4_L001_I1_001.fastq.gz
    ├── SAMPLENAME_S4_L001_R1_001.fastq.gz
    └── SAMPLENAME_S4_L001_R2_001.fastq.gz
input_fastqs_gex
└── SAMPLENAME
    ├── SAMPLENAME_S4_L001_I1_001.fastq.gz
    ├── SAMPLENAME_S4_L001_R1_001.fastq.gz
    └── SAMPLENAME_S4_L001_R2_001.fastq.gz
```

### Configure the pipeline

The pipeline must be appropriatly configured to your data. A detailed [Readme](config/README.md) can be found in the `config` directory. 


### Preprocessing

**IndexHopping removal**  
In case of combined GEX and ADT NovaSeq sequencing data, scripts are provided to clean up the data before a run. Please consult the [Readme](workflow/scripts/index_hopping_removal/README.md) here.


## Running the pipeline
Following the configuration of the pipeline a dryrun can be started using:
```
snakemake --use-conda --printshellcmds --dry-run
```


