# README testdata


## Example data

We provide [example data for a test run](https://drive.google.com/drive/folders/14clt2_E_P0-HEXlJwH1fHCk5KhpPpxMc?usp=share_link) with three hashed samples of human PBMC cells, so that hashing deconvolution, GEX analysis and ADT analysis can be performed. 

The raw test data comprises

- ADT FASTQ files
- GEX FASTQ files

All config files are ready to use in the `testdata` subdirectory.
- [HashingFile](HashingFile_PBMC_D1.csv) with hashtag barcodes
- [featureReferenceFile](feature_reference.txt) with all ADT barcodes 


### Quick test run

As an alternative we provide the count data that is generated on the example data described above, to allow skipping the resource-intensive cellranger count and CITE-Seq steps.

To start a quick test run (please refer to section `Installation instructions` in the [project readme](../README.md) for details on steps 1 and 2):

1) Install Snakemake, mamba and snakedeploy on your system
2) Deploy the workflow with the automated `snakedeploy` command that requires internet access, or using a local set up with cloning the repository, as described in the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules).
3) Unpack the test data matrices by running the following command in the gExcite working directory (usually `gExcite_pipeline`)
```
mv testdata/results_and_fastqs.tar.gz  . ; tar -xf results_and_fastqs.tar.gz
```

The directories `results` and `fastqs`, containing the raw count matrices, are now available in the working directory.

4) Submit a dry-run to test the configuration

```
snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --printshellcmds --dry-run --rerun-triggers mtime
```
NOTE: the parameter `--rerun-triggers mtime` makes sure only changes to the input data triggers a rerun of the pipeline.  

5) Start the Snakemake workflow

```
snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --printshellcmds --rerun-triggers mtime
```


### Full test run
To start a full test run that also includes the resource-intensive cellranger count and CITE-Seq steps (please refer to section `Installation instructions` in the [project readme](../README.md) for details on steps 3-5):

1) Download the FASTQ files archive and extract it with `unzip gexcite_testdata_fastqs.zip`
2) Move or link them into a subdirectory called `fastqs` in the gExcite working directory (usually `gExcite_pipeline`). Make sure you follow the [expected folder structure](../README.md) with one subdirectory per sample.
3) Install Snakemake, mamba and snakedeploy on your system.
4) Deploy the workflow with the automated `snakedeploy` command that requires internet access, or using a local set up with cloning the repository as described in the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules).
5) Install the [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) software. Follow the instructions on the 10xGenomics installation support page. Download the cellranger references as well.
6) Insert the paths to the available cellranger software and reference directory into the config file `config/config.yaml`
    - In section [`resources`], `reference_transcriptome` needs to point to the location of the genomic reference used for the cellranger mapping
    - In sections [`tools`][`cellranger_count_gex`] and [`tools`][`cellranger_count_adt`], `call` needs to point to the the path to the cellranger installation

Refer to the [config README file](../config/README.md) for more details

7) Submit a dry-run to test the configuration

```
snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --printshellcmds --dry-run --rerun-triggers mtime
```
NOTE: the parameter `--rerun-triggers mtime` makes sure only changes to the input data triggers a rerun of the pipeline.  

8) Start the Snakemake workflow with

```
snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --printshellcmds --rerun-triggers mtime
```
