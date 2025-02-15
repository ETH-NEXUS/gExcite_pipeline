# Configuration

Before using the pipeline the following files need to be provided/adapted:

1. config.yaml
2. samplemap
3. HashingFile
4. featureReferenceFile

## config.yaml

Before running the pipeline the `config.yaml` file needs to be adapted to contain the input and output paths for the intended analysis.

Adaptation necessary for a default run:

- In section [`inputOutput`] the input directories `input_fastqs_gex`, `input_fastqs_adt` need to point to the location of the respective FASTQ files. This location needs to be specified relative to the gExcite working directory (usually `gExcite_pipeline/`).
- In section [`resources`], `reference_transcriptome` needs to point to the location of the genomic reference used for the cellranger mapping
- In sections [`tools`][`cellranger_count_gex`] and [`tools`][`cellranger_count_adt`], `call` needs to point to the the path to the cellranger installation
- Section [`computingResources`] needs to list the resources that can be assigned to the analysis steps your data and batch system. The example resources specify memory per job, not thread.
- Section [`scampi`][`resources`] needs to be filled with the cell type information, selected genes to show in expression plots, and gene sets for the GSVA analysis.

The default [example config file](config.yaml) is pointing to the example data input files.

## samplemap

Further, a "samplemap" must be provided specifying sample-specific parameters in a tab-delimited text file.
A pre-configured samplemap ready to run on the test data that can be adapted is provided in this directory.

Example samplemap:

```
sample    HashingFile     nTargetCells    featureReferenceFile
sampleA   HashingFileA    10000           featureReferenceFileA
sampleB   HashingFileB    15000           featureReferenceFileB
```

With one line per set of samples

- `sample` contains the sample identifier that is used throughout the pipeline
- `HashingFile` contains the full path to the comma separated text file containing the hashtag barcodes and their assignment to individual sample names (see [HashingFile](#hashingfile)).

- `nTargetCells` corresponds to the number of targeted cells for the sample.
- `featureReferenceFile` corresponds to the ADT feature reference file for the sample set. For further information please consult the Cellranger tool documentation.  
**NOTE:** the path to the FeatureReferenceFile must be relative to the gExcite working directory (usually `gExcite_pipeline/`).

## HashingFile

In case of hashed samples, the hashtag barcodes, the hashtag names, and the corresponding sample names must be associated with the sample set. To do so, we need a comma-separated file with the following structure:

```
Barcode1,TagName1,sampleA
Barcode2,TagName2,sampleB
```

##  featureReferenceFile

The "featureReferenceFile" is a comma-separated text file describing all ADT antibodies used in the experiment at hand.
For further information please consult the Cellranger tool [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis#feature-ref).  
An example `feature_reference.txt` ready to run on the test data is available in the `testdata` directory.

