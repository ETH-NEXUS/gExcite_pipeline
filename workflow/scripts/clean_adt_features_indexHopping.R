# Clean index hopping output for ADT data
# singlecell_analysis
# Lourdes Rosano, January 2020

library(Matrix)

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

args = commandArgs(trailingOnly = TRUE)

features_in   = args[1] # features file from index hopping removal
matrix_in     = args[2] # matrix file from index hopping removal
adt_feat_in   = args[3] # adt reference features file, expects specific file structure
features_out  = args[4] # new (clean) features file
matrix_out    = args[5] # new (clean) matrix file

## Read in matrix file
# From: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices
mat = readMM(file = matrix_in)
print(dim(mat))
print(str(mat))

## Read in features file
feature_names = read.delim(features_in, header=FALSE, stringsAsFactors=FALSE, quote="", check.names=FALSE)
print(dim(feature_names))
print(str(feature_names))

# Assign features as rownames of the matrix
rownames(mat) = feature_names$V1
print(str(mat))

## Read in ADT reference features list
# Assumes csv file with the following header: id,name,read,pattern,sequence,feature_type
adt_features = read.csv(adt_feat_in, header=TRUE, stringsAsFactors=FALSE, quote="", check.names=FALSE)
print(dim(adt_features))
print(str(adt_features))
print(table(feature_names$V1 %in% adt_features$id))

new_features = subset(feature_names, subset = V1 %in% adt_features$id)
print(str(new_features))

# Sanity check that subsetting of ADT features yielded the expected result
stopifnot(nrow(new_features)==nrow(adt_features))
stopifnot(length(setdiff(new_features$V1,adt_features$id))==0)

## Subset matrix to only contain rows corresponding to ADT features (remove all GEX)
print(table(rownames(mat) %in% new_features$V1))

subMat = mat[rownames(mat) %in% new_features$V1,,drop=F]
print(dim(subMat))
print(str(subMat))

# Sanity check that subsetting of ADT matrix yielded the expected result
stopifnot(ncol(mat)==ncol(subMat))
stopifnot(nrow(subMat)==nrow(adt_features))
stopifnot(length(setdiff(adt_features$id,rownames(subMat)))==0)

rownames(subMat) = NULL
print(str(subMat))

## Correct type of feature for ADT data See: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices
# "The third column identifies the type of feature, which will be one of Gene Expression, Antibody Capture, CRISPR, or CUSTOM, depending on the feature type."
print(head(new_features))
print(table(new_features$V3))

new_features$V3 = "Antibody Capture"
print(head(new_features))
print(table(new_features$V3))

writeMM(subMat, file=matrix_out)
write.table(new_features, file=features_out, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
