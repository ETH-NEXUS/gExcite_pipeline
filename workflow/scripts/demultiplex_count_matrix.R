
# Script is filtering the cellranger runs based on the hashing barcode lists from the preprocessing pipeline.
# The output is written into the correct folder structure for the singlecell pipeline to get started.
library(optparse)
library(Matrix)
library(R.utils)

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

option_list <- list(
  make_option("--sampleMap", type = "character", help = "Path to the samplemap."),
  make_option("--sample", type = "character", help = "Sample Name."),
  make_option("--sampleSet", type = "character", help = "Sample Set Name.")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

sampleMap <- read.csv(opt$sampleMap, header = TRUE, sep = "\t")
print(sampleMap)

filter_cellRanger <- function(barcodeList, inputFolder, outputFolder, hashedSampleSet, sampleName) {
  # Read in hashing barcodes
  hashbarcodes <- read.table(barcodeList)
  # Read in cellranger barcodes
  sampleCellRanger <- read.table(paste(inputFolder, paste(hashedSampleSet, "barcodes.tsv", sep = "."), sep = "/"))
  print(hashbarcodes[1, 1])
  barcode_length <- nchar(toString(hashbarcodes[1, 1]))
  print(paste0("Barcode length:", barcode_length))

  # Convert cellranger barcodes into appropriate format if necessary(remove tag at the end)
  # See: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam
  # eg. AGAATGGTCTGCAT-1
  sampleCellRanger$V2 <- substr(sampleCellRanger$V1, 1L, barcode_length)
  print(head(sampleCellRanger))


  print(paste0("Detected cells by hashing: ", nrow(hashbarcodes)))
  print(head(hashbarcodes))
  print(paste0("Detected cells by CellRanger on pooledSample: ", nrow(sampleCellRanger)))
  print(head(sampleCellRanger))
  # Read in cellranger matrix
  mat <- Matrix::readMM(file = paste(inputFolder, paste(hashedSampleSet, "matrix.mtx", sep = "."), sep = "/"))
  print(dim(mat))
  print(str(mat))

  # Assign cellranger barcodes as colnames of the matrix
  colnames(mat) <- sampleCellRanger$V2
  print(str(mat))

  # Retrieve intersection of cell barcodes in index hopping and cellranger outputs
  intersection <- intersect(hashbarcodes$V1, sampleCellRanger$V2)
  # Sanity check that at least some barcodes could be intersected
  stopifnot(length(intersection) > 0)
  print(paste0("Final intersection cells: ", length(intersection)))
  # Remove excluded barcodes also from the matrix
  print(table(colnames(mat) %in% intersection))

  subMat <- mat[, colnames(mat) %in% intersection, drop = F]
  print(dim(subMat))
  print(str(subMat))

  # Sanity check that subsetting of ADT matrix yielded the expected result
  stopifnot(nrow(mat) == nrow(subMat))
  stopifnot(ncol(subMat) == length(intersection))
  stopifnot(length(setdiff(intersection, colnames(subMat))) == 0)

  # Write unzipped out files with samplename
  colnames(subMat) <- NULL
  print(str(subMat))
  outbarcodes <- paste(outputFolder, paste(hashedSampleSet, sampleName, "barcodes.tsv", sep = "."), sep = "/")
  print(outbarcodes)
  outmatrix <- paste(outputFolder, paste(hashedSampleSet, sampleName, "matrix.mtx", sep = "."), sep = "/")
  print(outmatrix)
  outfeatures <- paste(outputFolder, paste(hashedSampleSet, sampleName, "features.tsv", sep = "."), sep = "/")
  print(outfeatures)
  write.table(intersection, file = outbarcodes, row.names = FALSE, quote = FALSE, col.names = FALSE)
  Matrix::writeMM(subMat, file = outmatrix)
  file.copy(paste(inputFolder, paste(hashedSampleSet, "features.tsv", sep = "."), sep = "/"), outfeatures)
  # Write zipped outfiles without samplename mimiking  as input for Seurat
  zippFileDir <- paste(outputFolder, paste(hashedSampleSet,sampleName, sep="."), "outs/filtered_feature_bc_matrix/", sep = "/")
  dir.create(zippFileDir, recursive = TRUE)
  R.utils::gzip(outbarcodes, paste(zippFileDir, "barcodes.tsv.gz", sep = "/"), overwrite = TRUE, remove = FALSE)
  R.utils::gzip(outmatrix, paste(zippFileDir, "matrix.mtx.gz", sep = "/"), overwrite = TRUE, remove = FALSE)
  R.utils::gzip(outfeatures, paste(zippFileDir, "features.tsv.gz", sep = "/"), overwrite = TRUE, remove = FALSE)
}


hashedSampleSet <- opt$sampleSet
sample <- opt$sample
hashingFile = sampleMap$HashingFile[which(sampleMap$sample == opt$sampleSet)]
print(hashingFile)
sampleTagMap <- read.csv(hashingFile, header = FALSE, sep = ",")
colnames(sampleTagMap) <- c("barcode", "tagName", "sampleName")
# Create the required output directories
baseoutdir <- "results/"
outdirADT <- paste(baseoutdir, "cellranger_adt", sep = "/")
outdirGEX <- paste(baseoutdir, "cellranger_gex", sep = "/")
dir.create(outdirADT, showWarnings = FALSE)
dir.create(outdirGEX, showWarnings = FALSE)
# Read the existing input directory
hashingdir <- paste("results/pooled_samples/hashing_analysis/", sep = "/")
indirADT <- paste("results/pooled_samples/cellranger_adt/", sep = "/")
indirGEX <- paste("results/pooled_samples/cellranger_gex/", sep = "/")

tag  <- sampleTagMap$tagName[which(sampleTagMap$sampleName == sample)]

# Start defining input variables for CellRanger Filter Function
barcodeFile <- paste(hashedSampleSet, tag, sample, "barcodes_singlets.txt", sep = ".")
# Filter CellRanger
filter_cellRanger(paste(hashingdir, barcodeFile, sep = ""), indirADT, outdirADT, hashedSampleSet, sample)
filter_cellRanger(paste(hashingdir, barcodeFile, sep = ""), indirGEX, outdirGEX, hashedSampleSet, sample)

# Creating demultiplexing complete file
file.create(paste(outdirADT, "/demultiplexing_complete.",hashedSampleSet,".txt", sep=""))
file.create(paste(outdirGEX, "/demultiplexing_complete.",hashedSampleSet,".txt", sep=""))
