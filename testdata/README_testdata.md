### README testdat

For a quick test run starting after the resource-intensive Cellranger steps see next section "Quick test run".  



## Example data

We provide [example data for a test run](https://drive.google.com/drive/folders/14clt2_E_P0-HEXlJwH1fHCk5KhpPpxMc?usp=share_link) with three hashed samples of human PBMC cells, so that hashing deconvolution, GEX analysis and ADT analysis can be performed. 

The raw test data comprises

- ADT FASTQ files
- GEX FASTQ files

All config files are ready to use in the `testdata` subdirectory.
- HashingFile with hashtag barcodes
- featureReferenceFile with all ADT barcodes
- samplemap with sample information  


### Quick test run

A quick test run on the example data can be performed that starts after the resource-intensive cellranger count and CITE-Seq steps.  

To start a quick test run:

1) Install Snakemake, mamba and snakedeploy on your system
2) Deploy the workflow with the automated `snakedeploy` command that requires internet access, or using a local set up as described in the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules).
3) Install the software for the pipeline using Snakemake's `--use-conda` functionality
4) Unpack the test data matrices by running the following command in the gExcite working directory (usually `gExcite_pipeline`)
```
mv testdata/results_and_fastqs.tar.gz  . ; tar -xf results_and_fastqs.tar.gz
```

The directories `results` and `fastqs`, containing the raw count matrices, are now available in the working directory.

5) Start the Snakemake workflow

```
snakemake -s workflow/Snakefile_testdata --configfile testdata/config_testdata.yaml --use-conda --printshellcmds
```


### Full test run
To start a full test run that also includes the resource-intensive cellranger count and CITE-Seq steps:

1) Download the FASTQ files
2) Move or link them into a subdirectory called `fastqs` in the gExcite working directory (usually `gExcite_pipeline`). Make sure you follow the [expected folder structure](README.md) with a subdirectory per sample.
3) Install Snakemake, mamba and snakedeploy on your system
4) Deploy the workflow with the automated `snakedeploy` command that requires internet access, or using a local set up as described in the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules).
5) Install the software for the pipeline using Snakemake's `--use-conda` functionality
6) Install the [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)software. Follow the instructions on the 10xGenomics installation support page. Download the cellranger references as well.
4) Insert the paths to the available cellranger software and reference directory into the testdata config `testdata/config_testdata.yaml`
    - `cellranger_count_gex`
    - `cellranger_count_adt`
    - `reference_transcriptome`

4. Start the Snakemake workflow with

```
snakemake -s workflow/Snakefile_testdata --configfile testdata/config_testdata.yaml --use-conda --printshellcmds
```
