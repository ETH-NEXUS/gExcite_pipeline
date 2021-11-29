# Clean detected cells from index hopping output for ADT data
# singlecell_analysis
# Lourdes Rosano, January 2020

library(Matrix)

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

args = commandArgs(trailingOnly = TRUE)

barcodes_in     = args[1] # barcodes file from index hopping removal
matrix_in       = args[2] # matrix file from index hopping removal
cr_barcodes     = args[3] # cellranger filtered barcodes file
barcode_length  = args[4] # barcode length. default = 16
barcodes_out    = args[5] # new (clean) barcodes file
matrix_out      = args[6] # new (clean) matrix file

# Read in index hopping barcodes
sampleIndex = read.table(barcodes_in)
# Read in cellranger filtered barcodes
sampleCellRanger = read.table(cr_barcodes)

print(paste0("Detected cells by index hopping removal method: ", nrow(sampleIndex)))
print(head(sampleIndex))
print(paste0("Detected cells by CellRanger: ", nrow(sampleCellRanger)))
print(head(sampleCellRanger))

# Read in index hopping matrix
mat = readMM(file = matrix_in)
print(dim(mat))
print(str(mat))

# Assign index hopping barcodes as colnames of the matrix
colnames(mat) = sampleIndex$V1
print(str(mat))

barcode_length = as.integer(barcode_length)
print(paste0("Barcode length:", barcode_length))

# Convert cellranger barcodes into appropriate format (remove tag at the end)
# See: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam
# eg. AGAATGGTCTGCAT-1
sampleCellRanger$V2 <- substr(sampleCellRanger$V1, 1L, barcode_length)
print(head(sampleCellRanger))

# Retrieve intersection of cell barcodes in index hopping and cellranger outputs
intersection <- intersect(sampleIndex$V1,sampleCellRanger$V2)
# Sanity check that at least some barcodes could be intersected
stopifnot(length(intersection) > 0)
print(paste0("Final intersection cells: ", length(intersection)))

# Remove excluded barcodes also from the matrix
print(table(colnames(mat) %in% intersection))

subMat = mat[,colnames(mat) %in% intersection,drop=F]
print(dim(subMat))
print(str(subMat))

# Sanity check that subsetting of ADT matrix yielded the expected result
stopifnot(nrow(mat)==nrow(subMat))
stopifnot(ncol(subMat)==length(intersection))
stopifnot(length(setdiff(intersection,colnames(subMat)))==0)

colnames(subMat) = NULL
print(str(subMat))

write.table(intersection, file = barcodes_out, row.names=FALSE, quote=FALSE, col.names=FALSE)
writeMM(subMat, file=matrix_out)
