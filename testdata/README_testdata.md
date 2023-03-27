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

As an alternative we provide the count data that is generated on the example data described above, to allow skipping the resource-intensive cellranger count and CITE-Seq steps. Note that the test run requires approximately 12 GB of memory to complete successfully.

To start a quick test run

1) Clone the gExcite git repository and go into the new directory `gExcite_pipeline` that will be referred to as "the gExcite working directory" in this documentation.

```
git clone https://github.com/ETH-NEXUS/gExcite_pipeline.git ;
cd gExcite_pipeline
```

2) Unpack the test data matrices by running the test data preparation bash script in the gExcite working directory.

```
sh prepare_quick_testrun.sh
```

The directories `results` and `fastqs`, containing the raw count matrices, are now available in the working directory.

3) Install Snakemake and mamba on your system

```
conda create -c bioconda -c conda-forge --name snakemake mamba snakemake ;
conda activate snakemake
```

4) Do a dry-run to test the configuration

```
snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --printshellcmds --dry-run --rerun-triggers mtime
```

NOTE: the parameter `--rerun-triggers mtime` makes sure only changes to the input data trigger a rerun of the pipeline.  

5) Start the Snakemake workflow

```
snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --printshellcmds --rerun-triggers mtime --cores 1
```

**NOTE:** if the pipeline should be run on a compute cluster using a job scheduling system (e.g. LSF, Slurm) the command needs to be adjusted accordingly. Please refer to the [Snakemake documentation on cluster execution](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) for platform-specific details.

### Full test run

To start a full test run that also includes the resource-intensive cellranger count and CITE-Seq steps:

1) Clone the gExcite git repository and go into the new directory `gExcite_pipeline` that will be referred to as "the gExcite working directory" in this documentation.

```
git clone https://github.com/ETH-NEXUS/gExcite_pipeline.git ;
cd gExcite_pipeline
```

2) Download the FASTQ files archive and extract it with `unzip gexcite_testdata_fastqs.zip`

3) Move the extracted directory (named "fastqs") into the gExcite working directory. 

4) Install Snakemake and mamba on your system

```
conda create -c bioconda -c conda-forge --name snakemake mamba snakemake ;
conda activate snakemake
```

5) Install the [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) software. Follow the instructions on the 10xGenomics installation support page to install cellranger and to include the cellranger binary to your path. Download the cellranger references.
6) Insert the paths to the available cellranger software and reference directory into the config file `config/config.yaml`
    - In section [`resources`], [`reference_transcriptome`] needs to point to the location of the genomic reference used for the cellranger mapping
    - In sections [`tools`][`cellranger_count_gex`] and [`tools`][`cellranger_count_adt`]:  
    [`call`] needs to point to the path to the cellranger installation  

7) Do a dry-run to test the configuration

```
snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --printshellcmds --dry-run
```

8) Start the Snakemake workflow

```
snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --printshellcmds --cores 1
```

**NOTE:** if the pipeline should be run on a compute cluster using a job scheduling system (e.g. LSF, Slurm) the command needs to be adjusted accordingly. Please refer to the [Snakemake documentation on cluster execution](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) for platform-specific details.
