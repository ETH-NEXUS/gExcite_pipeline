[![Snakemake](https://img.shields.io/badge/snakemake-≥6.12.1-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![Code style: snakefmt](https://img.shields.io/badge/code%20style-snakefmt-000000.svg)](https://github.com/snakemake/snakefmt)
[![GitHub Super-Linter](https://github.com/ETH-NEXUS/scGATE_workflow/workflows/Lint%20Code%20Base/badge.svg)](https://github.com/marketplace/actions/super-linter)
# scGate 

## General overview


## Installation instructions
### Pipeline
Given conda is installed on you system the pipeline can be set up using snakedeploy

First create and activate an environment including Mamba, Snakemake and Snakedeploy:

```
conda create -c bioconda -c conda-forge --name snakemake mamba snakemake snakedeploy
conda activate snakemake
```

Snakedeploy can now be used to deploy the workflow. 

```
snakedeploy deploy-workflow https://github.com/ETH-NEXUS/scGATE_workflow --tag master .
```

### Dependencies
Most of the software used in the default workflow can be installed in an automated fashion using snakemake's [--use-conda](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management) functionality. 
The following software needs to be installed manually.

- [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger): Follow the instructions on the 10xGenomics installation support page to install cellranger and to include the cellranger binary to your path.
Webpage: [https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation)

- [CITE-seq-Count](https://hoohm.github.io/CITE-seq-Count/) 
 
```
pip install CITE-seq-Count
```


## Example data



## Before running the pipeline
Before the pipeline can run the ADT and the GEX fastq files have to be placed in a specific folder structure, the pipeline has to be configured appropriately and a preprocessing has to be performed if necessary. 

### Prepare the fastq files
The pipeline expects the fastq files per sample to be in the following folder structure and adhere to the naming schema: 

/path/to/**input_fastqs_gex_directory**/**SAMPLENAME**/**SAMPLENAME**_S[Number]_L00[Lane Number]_[Read Type]_001.fastq.gz


### Configure the pipeline

The pipeline must be appropriatly configured to your data. A detailed [Readme](config/README.md) can be found in the 'config' directory. 


### Preprocessing

**IndexHopping removal**  
In case of combined GEX and ADT NovaSeq data scripts are provided to clean up the data before a run. Please consult the [Readme](workflow/scripts/index_hopping_removal/README.md) here.


## Running the pipeline
Following the configuration of the pipeline a dryrun can be started using:
```
snakemake --use-conda --printshellcmds --dry-run
```


