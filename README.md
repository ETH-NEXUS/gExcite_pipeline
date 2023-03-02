[![Snakemake](https://img.shields.io/badge/snakemake-≥6.12.1-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

# gExcite pipeline

![Workflow Figure](https://github.com/ETH-NEXUS/gExcite_pipeline/blob/main/images/Workflow_Figure_gExcite.png)

## General overview

gExcite is a start-to-end workflow embedded in Snakemake that provides both, gene expression and CITE-seq analysis, as well as hashing deconvolution.  
For an overview of all steps please see the Snakemake [rulegraph](https://github.com/ETH-NEXUS/gExcite_pipeline/blob/main/images/gExcite_pipeline_rulegraph.png).

## Remark

This workflow makes use of Snakemake's functionality to include external workflows as a [module](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#snakefiles-modules).
scAmpi, a workflow that provides basic scRNA processing steps, is included as a module into gExcite. Note that all documentation regarding scAmpi (especially regarding config file entries that must be adapted depending on the disease) can only be found in the ![scAmpi](https://github.com/ETH-NEXUS/scAmpi_single_cell_RNA) git repository.

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

Note: Snakemake needs to access the internet for this set up. With Snakemake 7.13 there is also support for a local set up of modules. Please refer to the [Snakemake documentation on modules](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules) for more details.

### Dependencies

Most of the software used in the default workflow can be installed in an automated fashion using Snakemake's [--use-conda](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management) functionality.
The following software needs to be installed manually.

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

### Configure the pipeline

The pipeline must be appropriately configured to your data. A detailed [Readme](config/README.md) can be found in the `config` directory.

### Optional preprocessing

**IndexHopping removal**  
In case of combined GEX and ADT NovaSeq sequencing data, scripts are provided to clean up the data before a run. Please consult the [Readme](workflow/scripts/index_hopping_removal/README.md) here.

## Running the pipeline

Following the configuration of the pipeline a run can be started using:

```
# dry run
snakemake --use-conda --printshellcmds --dry-run
# analysis run
snakemake --use-conda --printshellcmds
```

## Example data

We provide [example data for a test run](https://drive.google.com/drive/folders/14clt2_E_P0-HEXlJwH1fHCk5KhpPpxMc?usp=share_link) with three hashed samples of human PBMC cells. With this data hashing deconvolution, GEX analysis and ADT analysis can be performed.  
The test data comprises

- ADT FASTQ files
- GEX FASTQ files

HashingFile with hashtag barcodes, featureReferenceFile with all ADT barcodes, and samplemap are available in the `testdata` directory.  

To start a test run

1) Download the FASTQ files
2) Move or link them into a subdirectory called `fastqs` in the gExcite working directory (usually `gExcite_pipeline`)
3) Follow the software installation instructions
4) Insert the paths to the available cellranger software and reference transcriptome into the testdata config `testdata/config_testdata.yaml`
    - `cellranger_count_gex`
    - `cellranger_count_adt`
    - `reference_transcriptome`

4. Start the Snakemake workflow with

```
snakemake -s workflow/Snakefile_testdata --configfile testdata/config_testdata.yaml --use-conda --printshellcmds
```

## Quick test run

A quick test run on the example data can be performed that starts after the resource-intensive cellranger count and CITE-Seq steps.  

To start a quick test run:

1) Follow the software installation instructions
2) To unpack the test data matrices run the following command in the gExcite working directory (usually gExcite_pipeline)

```
mv testdata/results_and_fastqs.tar.gz  . ; tar -xf results_and_fastqs.tar.gz
```

The directories `results` and `fastqs`, containing the raw count matrices, are now available in the working directory.

3) Start the Snakemake workflow

```
snakemake -s workflow/Snakefile_testdata --configfile testdata/config_testdata.yaml --use-conda --printshellcmds
```
