# Overview
The folder contains scripts dedicated to the aggregation, normalization and differential expression analysis of single cell data.

All necessary dependencies are available within the conda environment described in `envs/de_analysis.yml`.
In order to use the environment please use the yaml file to create the environment and then activate it.
- Install conda on your machine (https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- create the conda environment: `conda env create -f de_analysis.yml`
- activate the environment: `conda activate de_analysis`

All scripts have a dedicated detailed documentation at the beginning of the files, that explains functionality, input, output, dependencies and other details. It is strongly suggested to familiarize with the specific documentation before running the scripts for the first time.

When working with Single Cell data, we provide two separate methods for differential expression analysis: DESeq2 or Seurat findMarkers.

DESeq2 framework is based on single sample data and works by "bulkyfying" the single-cell data available. In a nutshell, counts are aggregated, so that for each feature (e.g. cell type), all single-cell counts are merged in bulk, per sample counting. When this is done, the assumption is that the data behave similarly to bulk-rna data and can be processed with standard DESeq procedures. This method has the following advantages:
The advantages are:
- flexibility of the DE contrasts, as the user can potentially apply any complex contrast allowed by the DESeq2 workflow and the script allows for simplified definition of contrasts with multiple arbitrary covariates
- using a common, widely known and accepted tool for DE
- the aggregation and DE are performed starting from the single samples and not from an integrated Seurat object. This ensures
The disadvantages are:
- low quality results on ADT data

The Seurat R package provides a function called findMarkers that is built to perform differential expression analysis on single-cell data specifically. The method is mature enough to be suggested in Seurat's vignette for both GEX and ADT data.
The advantages are:
- working on a single, aggregated seurat object containing all necessary results simplifies the workflow and allows for easier manipulation and subsetting of the data
- better overall result quality on ADT data
The disadvantages are:
- The model allows only pairwise comparisons between feature values without covariates. Complex contrasts are not possible at the moment

# Installation
All scripts are programmed to be run through the command line instead of in the R interpreter. All script should already have "execute" permissions and can be run by simply calling the script itself.
Options must be provided through the command line. A list of options with detailed explanation is available in each script's documentation. A quick reference can be retrieved by using the option --help

The sripts have been tested on the following software and packages versions:
dbscan 1.1_10
DESeq2 1.32.0
dynamicTreeCut 1.63-1
ggplot2 3.3.5
ggrepel 0.9.1
GSVA 1.40.1
multcomp 1.4-17
optparse 1.6.6
pheatmap 1.0.12
R 4.1.0
RColorBrewer 1.1-2
reshape2 1.4.4
scran 1.20.1
Seurat 4.0.4
SingleCellExperiment 1.14.1
UpSetR 1.4.0
uwot 0.1.10
zeallot 0.1.0

# Preprocessing
The differential expression procedure depends on the results of the GEX analysis. For each sample, a dedicated RDS file containing a SingleCellExperiment object is necesary.
The position of the RDS files and the name to be used during the analysis must be provided by the user through a metadata table. The metadata table is a tsv table with two column and no header. The first column must contain the path to the RDS files to aggregate, one file for each row. The second column must contain the sampleID to use within the analysis for each specific file. Following is an example of metadata file

/path/to/sample1/sample1.RDS    sample1ID
/path/to/sample2/sample2.RDS    sample2ID

DESeq2-based differential expression analysis requires the user to define the contrast of interest. This information must be provided in a text file where all variables of interest are provided in the same line, separated by tabulations. The script will treat the last variable as the main contrast variable and any previous variable as covariates. All variable names must match the variable names of the provided SingleCellExperiment objects. Following is an example of contrast file:

stage   biopsy_localisation

The standalone normalization script requires clincal data for the samples. Clinical data must be provided as a tab-delimited table with header, one sample per line. The table must contain a column storing the names of each sample. These sample names must match the sample names available in the aggregated RDS files (column names)

#cohort_aggregation.R
The aggregation script is built to load a set of standard RDS files output from the scPipeline. Each RDS file contains data for a single sample and all samples used in a run of this script should come from the same cohort.
The script completes the aggregation by, in fact, aggregating all samples by the grouping variable of interest. Examples of grouping variables can be "celltype_major" or "ct_cl"

#cohort_group_normalize_cluster_plot.R
The normalization script is built to load the results from the aggregation step and it assumes the output structure and filenames used by the aggregation script.
The script loads the aggregated data and performs: variance stabilization, rlog normalization and multiple QC plots.
Please note that this step IS NOT REQUIRED when doing a DE analysis, because THE DE SCRIPT ALREADY HAS A BUILT-IN NORMALIZATION STEP. Normalizing twice is not suggested.

#cohort_group_deseq_plot.R
The differential expression script is built to load the results from the aggregation step and it assumes the output structure and filenames used by the aggregation script.
The script loads the aggregated data and performs a differential expression analysis using DESeq2. Contrasts are based on a user-provided clinical data table. The script provides the user with the DE table, plots about the gene expression and plot about significant genes (if any). Missing plots for a specific condition imply no significant gene was found

#de_findmarkers.R
Alternative differential expression script that makes use of n integrated Seurat object containing counts and metadata for the project. This script is to be preferred over DESeq2 when working with ADT or zero-inflated data, but has the limitation to be able to work only with pairwise contrasts without covariates. Please refer to the document "readme_differential_expression" in the folder "doc" for a more in-depth explanation on the pros and cons of each method

#findmarkers_plotheatmap.R
Script created to be attached downstream of the DE analysis. Standalone in order not to clutter the option list.
It takes in input the DE output from DE_findMarkers.R and creates an average-expression heatmap of the top ranked markers.

# Postprocessing
No post-processing is necessary for the results output by the differential expression analysis
