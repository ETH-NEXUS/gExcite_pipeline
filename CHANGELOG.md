# Changelog


## [1.0.5] - 2024-10-09 
- removed `SeqRunName` variable, as it is not needed by cellranger anymore

## [1.0.4] - 2023-09-19

### Changed
- with `Snakefile_no_hashing.smk` have the option to run gExcite without hashing deconvolution

### Fixed
- catch special case in script `analyse_citeseq.R` where cell barcodes have suffix "-1"

## [1.0.3] - 2023-08-21

### Changed
- Give cellranger ADT and cellranger GEX rule the number of threads specified in the config file (`--localcores`).
- in script `analyse_citeseq.R` have new parameter `--number_pca_adt`. With this parameter in the config file the number of PCA dimensions used for the UMAP calculation based on ADT counts can be adjusted. The number of PCA dimensions cannot be larger than the number of ADTs in the experiment at hand, otherwise the script fails.

## [1.0.2] - 2023-03-22
- Update version of scAmpi that is used as a module. scAmpi v2.0.7 contains fixes regarding R package discrepancies that caused various scripts to fail (using `ggsave` function).

## [1.0.1] - 2023-03-10
- Update test run settings to have all config files available and ready to use
- expand test run settings with the possibility to omit the resource-intensive cellranger steps
- update README.md and have separate testdata/README_testdata.md for describing the test run steps
- adapt the FASTQ input interpretation: specifying the input relative to the working directory is now mandatory, not merely recommended
- Improved documentation on the use of scAmpi as a module in gExcite

## [1.0.0] - 2022-12-08

First full, publicly available version of the gExcite analysis pipeline.


