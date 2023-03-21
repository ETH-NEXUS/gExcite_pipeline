[![Snakemake](https://img.shields.io/badge/snakemake-≥6.12.1-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

# gExcite pipeline

![Workflow Figure](https://github.com/ETH-NEXUS/gExcite_pipeline/blob/main/images/Workflow_Figure_gExcite.png)

## General overview

gExcite is a start-to-end workflow embedded in Snakemake that provides both, gene expression and CITE-seq analysis, as well as hashing deconvolution. The workflow is compatible with and tested on Linux only, other Unix systems (including MacOS) are currently not supported.
For an overview of all steps please see the Snakemake [rulegraph](https://github.com/ETH-NEXUS/gExcite_pipeline/blob/main/images/gExcite_pipeline_rulegraph.png).

## Remark

This workflow makes use of Snakemake's functionality to include external workflows as a [module](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#snakefiles-modules).
scAmpi, a workflow that provides basic scRNA processing steps, is included as a module into gExcite. Note that all documentation regarding scAmpi (especially regarding config file entries that must be adapted depending on the disease) can only be found in the [scAmpi](https://github.com/ETH-NEXUS/scAmpi_single_cell_RNA) git repository.

## Example data

We provide [example data for a test run](https://drive.google.com/drive/folders/14clt2_E_P0-HEXlJwH1fHCk5KhpPpxMc?usp=share_link) with three hashed samples of human PBMC cells, so that hashing deconvolution, GEX analysis and ADT analysis can be performed.
For more details see the [README](testdata/README_testdata.md) in the testdata subdirectory.

### Quick test run

A quick test run on the example data can be performed that starts after the resource-intensive cellranger count and CITE-Seq steps.  
For more details see the [README](testdata/README_testdata.md) in the testdata subdirectory.

## Installation instructions

### Pipeline

Given conda is installed on your system the pipeline can be set up using `snakedeploy`.

First, create and activate an environment including Mamba, Snakemake and Snakedeploy:

```
conda create -c bioconda -c conda-forge --name snakemake mamba snakemake snakedeploy ;
conda activate snakemake
```

Snakedeploy can now be used to deploy the workflow:

```
snakedeploy deploy-workflow https://github.com/ETH-NEXUS/gExcite_pipeline --tag main .
```

Note: Snakemake needs to access the internet for this set up. With Snakemake 7.13 there is also support for a local set up of modules. Please refer to the [Snakemake documentation on modules](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules) for more details.

### Dependencies

Most of the software used in the default workflow can be installed in an automated fashion using Snakemake's [--use-conda](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management) functionality when running the pipeline.
In case you would like to start from the raw sequencing data using cellranger processing, the following software needs to be installed manually.

- [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger): Follow the instructions on the 10xGenomics installation support page to install cellranger and to include the cellranger binary to your path.
Webpage: [https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation)

## Before running the pipeline

Before the pipeline can be run make sure that

- ADT and the GEX FASTQ files are provided in the folder structure specified below
- the pipeline config.yaml is configured to your data
- optional preprocessing was performed if necessary

### Prepare FASTQ files

The pipeline expects the FASTQ files per sample to be in the following folder structure, adhering to the naming schema:

/path/to/**input_fastqs_gex**/**SAMPLENAME**/**SAMPLENAME**\_S[Number]\_L00[Lane Number]\_[Read Type]\_001.fastq.gz  

Example:

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

Starting processing after the resource-intensive cellranger count and CITE-Seq steps requires the presence of dummy files in place of the FASTQ files described in this chapter. For more details see the [README](testdata/README_testdata.md) in the testdata subdirectory. For a full example of the required folder structure refer to the `results_and_fastqs.tar.gz` in the testdata subdirectory.

### Configure the pipeline

The pipeline must be appropriately configured to your data. A detailed [README](config/README.md) can be found in the `config` directory.

### Optional preprocessing

**IndexHopping removal**  
In case of combined GEX and ADT NovaSeq sequencing data, scripts are provided to clean up the data before a run. Please consult the [README](workflow/scripts/index_hopping_removal/README.md) here.

## Running the pipeline

Following the configuration of the pipeline a run can be started using:

```
# dry run
snakemake --use-conda --printshellcmds --dry-run
# analysis run
snakemake --use-conda --printshellcmds --cores 1
```
