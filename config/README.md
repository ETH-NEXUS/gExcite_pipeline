# Configuration

Before using the pipeline the following files need to be provided/adapted:

1. config.yaml
2. samplemap
3. HashingFile
4. featureReferenceFile

## config.yaml

Before running the pipeline the `config.yaml` file needs to be adapted to contain the input and output paths for the intended analysis.  
Adapt in

- the first section [`inputOutput`] the input directories `input_fastqs_gex`, `input_fastqs_adt`
- [`resources`][`reference_transcriptome`] the genomic reference used for the cellranger mapping
- [`tools`][`cellranger_count_gex`] and [`cellranger_count_adt`] the path to the cellranger installation
- [`computingResources`] the resources to your data and batch system. The example resources specify memory per job, not thread.
- [`scampi`][`resources`] the cell type information, selected genes to show in expression plots, and gene sets for the GSVA analysis.

## samplemap

Further, a "samplemap" must be provided specifying sample-specific parameters in a tab-delimited text file.
An example file ready for adaptation is provided in this directory.

Example:

```
sample    HashingFile     SeqRunName      nTargetCells    featureReferenceFile
sampleA   HashingFileA    SeqRunNameA     10000           featureReferenceFileA
sampleB   HashingFileB    SeqRunNameB     15000           featureReferenceFileB
```

- HashingStatus_x corresponds to a File containing the necessary information (see [HashingFile](#hashingfile)).
- SeqName_ADT_x corresponds to the sequencing sample name of the ADT sample for set 'x', as this parameter is only required for the CellRanger run of ADT data. It can be retrieved from the fastqs file names as follows:

```
[SequencingName_ADT_x]_S[Number]_L00[Lane Number]_[Read Type]_001.fastq.gz
```

Where Read Type is one of: I1, R1, R2.

- nTargetCells corresponds to the number of targeted cells for the sample. This parameter usually has to be specified for hashing experiments.
- featureReferenceFile corresponds to the ADT feature reference file for sample set 'x'. For further information please consult the cellranger tool documentation.

When parameters in the third and fourth column do not need to be provided, "." can be used instead.

## HashingFile

In case of hashed samples, the hashtag barcodes, the hashtag names, and the corresponding sample names must be associated with the sample set. To do so, we need a comma-separated file with the following structure:

```
Barcode1,TagName1,SampleA
Barcode2,TagName2,SampleB
```

##  featureReferenceFile

The "featureReferenceFile" is a comma-separated text file describing all ADT antibodies used in the experiment at hand.
For further information please consult the cellranger tool [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis#feature-ref).
