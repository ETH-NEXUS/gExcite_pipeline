# scGate 

## General overview


## Installation instructions
Most of the software used in the default workflow can be installed in an automated fashion using snakemake's [--use-conda](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management) functionality. 
The following software needs to be installed manually.

- [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger): Follow the instructions on the 10xGenomics installation support page to install cellranger and to include the cellranger binary to your path.
Webpage: [https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation)

- [CITE-seq-Count](https://hoohm.github.io/CITE-seq-Count/) 
 
```
>> pip install CITE-seq-Count
```


## Example data



## Before running the pipeline

#### Configure the pipeline

The pipeline must be appropriatly configured to your data. A detailed [Readme](config/README.md) can be found in the 'config' directory. 


#### Preprocessing

**IndexHopping removal**  
In case of combined GEX and ADT NovaSeq data scripts are provided to clean up the data before a run. Please consult the [Readme](workflow/scripts/index_hopping_removal/README.md) here.


## Running the pipeline
Following the configuration of the pipeline a dryrun can be started using:
```
snakemake --use-conda --printshellcmds --dry-run
```


