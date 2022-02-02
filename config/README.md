# Configuration

Before using the pipeline the following files need to be provided. 
1. config
2. samplemap
3. HashingFile
4. feature Reference File
 

## Configuration File
Before running the pipeline the `config` file needs to be adapted to contain the input and output paths for the intended analysis. Those are provided in the first section (`inputOutput`) of the config file. In addition to input and output paths, the reference must point to the genomic reference used for the cellranger mapping. Finally the correct call to cellranger and CITE-seq-Count needs to be included.

## SampleMap
Further, a "sample_map" must be provided to specify sample-specific parameters.
Am example file ready for adpatation is provided in this directory.

Example:
```
sample  HashingFile   SeqRunName      nTargetCells    featureReferenceFile
sampleA HashingFileA    SeqRunNameA     10000   featureReferenceFileA
sampleB HashingFileB    SeqRunNameB     15000   featureReferenceFileB
```

- HashingStatus_x corresponds to a File containing the necessary information (see [Hashing Info File](#hashing-info-file)).
- SeqName_ADT_x corresponds to the sequencing sample name of the ADT sample for set 'x', as this parameter is only required for the CellRanger run of ADT data. It can be retrieved from the fastqs file names as follows:
```
[SequencingName_ADT_x]_S[Number]_L00[Lane Number]_[Read Type]_001.fastq.gz
```
Where Read Type is one of: I1, R1, R2.
- nTargetCells_x corresponds to the original number of targeted cells for the sample in the given experiment. This parameter usually has to be specified for hashing experiments.
- featFile_ADT_x corresponds to the ADT feature reference file for sample set 'x'. For further information please consult the cellranger tool documentation.

When parameters in the third and fourth column do not need to be provided, "." can be used instead.

## Hashing Info File
In case of Hashed Samples we need to associate the hashtag barcodes, the tag names and the corresponting samplenames with the sampleset. To do so we need a file having the following structure:

```
Barcode1,TagName1,SampleA
Barcode2,TagName2,SampleB
```

## Feature Reference File

For further information please consult the cellranger tool [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis#feature-ref).
