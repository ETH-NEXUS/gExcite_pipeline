# Overview
To perform index hopping removal a rough guideline is provided here. 

## Preprocessing
First cellranger needs to be run on all samples that were sequenced together as if the snakemake workflow would run on these files. Note that ADT and GEX data needs to have seperate cellranger runs in order for this approach to work.  

## clean_and_filter_swapped.R
Needs the following input:
1. Path to the out Rdata file.
2. Path to a folder where all outfiles can be stored.
3. FDR threshold for cell calling (e.g. 0.01)
4. library size threshold for cell calling (e.g. 500)

## ADT Data
Unfortunately the output format of the above script is not a valid format to be processed with the pipeline therefore the next two scripts have to be run to get an appropriatly formated matrix. The order of the input files can be seen in the script. 
### clean_adt_barcodes_indexHopping.R
Output is a clean barcode and a clean matrix file. 
### clean_adt_features_indexHopping.R
Output is a clean feature file. 

## GEX Data
Unfortunately the ouput feature file of the above script contains Ensembl Ids instead of Gene Symbols. Using the following script helps to regenerate the appropriately formated feature file. 
### map_geneSymbols_indexHopping.R 


## Postprocessing
After above files have been generated, the pipeline can be started by placing the corrected matrix file into the appropriate cellranger folder of the pooled samples. Snakemake should then recognise these files as already generated and start without creating an uncorrected cellranger matrix. 
