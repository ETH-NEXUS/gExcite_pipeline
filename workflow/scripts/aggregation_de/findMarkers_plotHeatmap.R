#!/usr/bin/env Rscript
#################################################
## File name: DE_findMarkers.R
## Author: Matteo Carrara
## Date created: Oct 2021
## R Version: 4.1.0
##################################################

## TESTED ON THE FOLLOWING VERSIONS OF THE DEPENDENCIES:
## Seurat 4.1.0
## ggplot2 3.3.5
## ggrepel 0.9.1
## RColorBrewer 1.1-2
## pheatmap 1.0.12
## optparse 1.7.1

## GENERAL:
## This script creates an heatmap based on a set of DE results coming from the script DE_findMarkers.R. Due to the level of complexity, the script is standalone to allow more granular options.
## The heatmap that is build is based on a single specific contrast of interest. From that contrast, the script extracts all available RDS objects containing the DE results for all splitFeature values available and plots an heatmap of the average expression of the top N markers in all DE results.
## We can use as an example a DE analysis in which we split the data by celltype and used condition as contrast. We pass the options to the script to work on the contrast "treated vs untreated". The script will load all DE results for the aforementioned contrast for all celltypes. For each celltype it will extract the top N markers (based on the adjusted pvalue). The union of all markers selected in this way is used for the heatmap, by extracting the average expressions.

## ADDITIONAL IMPORTANT INFORMATION:
## This script takes as input a set of RDS objects generated by the DE analysis performed with the script DE_findMarkers.R and assumes all DE RDS objects are in the provided folder
## All options are required and must have a value set. In the future there might be a migration away from optparse and to a package allowing to
## automatically check if a required option has a value or not.

## INPUT:
## - The path where to find the RDS file storing the Seurat object containing the count data. It must be the same Seurat file provided for the DE analysis
## - The path to the folder that contains all input DE RDS objects for the contrast of interest
## - The name of the contrast of interest. It must be provided as seen in the filename of the output of DE_findMarkers.R. The filename structure is always <splitFeature>_<assay>-<contrast>
## - The name of the assay of interest. It must be provided as seen in the filename of the output of DE_findMarkers.R. The filename structure is always <splitFeature>_<assay>-<contrast>
## - The feature used to compute contrasts. It must correspond to the same option provided to DE script
## - The feature used to split the Seurat object. It must correspond to the same option provided to the DE script
## - The path to the folder where to save the output
## - Whether the count data in the main seurat object is log-normed or not. This information is part of the processing of the Seurat object and is very important to retrieve the correct values to show

## OUTPUT:
## A png image containing a heatmap. The heatmap shows the counts for the union of all top markers selected during DE for the contrast of interest in all the feature values used to split the Seurat object

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(Seurat)
  library(RColorBrewer)
  library(pheatmap)
  library(optparse)
})

option_list <- list(
  make_option("--seuratFile", type = "character", help = "Path to the RDS file containing the Seurat object with the counts and metadata"),
  make_option("--inputFolder", type = "character", help = "Path to the folder containing the all DE RDS files for the contrast of interest"),
  make_option("--contrast", type = "character", help = "Name of the contrast of interest as shown in the RDS objects' filenames"),
  make_option("--assay", type = "character", help = "Name of the assay of interest as shown in the RDS objects' filenames"),
  make_option("--contrastFeature", type = "character", help = "Name of the feature to use to compute the contrasts. It must be an existing metadata column in the Seurat object"),
  make_option("--splitFeature", type = "character", default = "", help = "Name of the feature to use to split the Seurat object. It must be an existing metadata column in the Seurat object"),
  make_option("--topMarkersNum", type = "integer", help = "Number of top markers to choose from each DE result table"),
  make_option("--outdir", type = "character", "Path to output directory"),
  make_option("--isLogNorm", type = "logical", default = FALSE, help = "Is the data stored in the Seurat object log-normalized?"),
  make_option("--pvalueadjThreshold", type = "double", default = NULL, help = "Is the selection of features to plot done on all genes (NULL) or only on genes under a certain PValue adjusted threshold?")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

"%&%" <- function(a, b) paste(a, b, sep = "")

seuratFile <- opt$seuratFile
inputFolder <- opt$inputFolder
mycontrast <- opt$contrast
splitFeature <- opt$splitFeature
contrastFeature <- opt$contrastFeature
myassay <- opt$assay
mytop <- opt$topMarkersNum
outdir <- opt$outdir %&% "/"
isLogNorm <- opt$isLogNorm
mypval <- opt$pvalueadjThreshold

message("Loading the Seurat object")
seurat_subset <- readRDS(seuratFile)

message("Loading the DE results")
myfiles <- list.files(inputFolder)
myline <- myassay %&% "_" %&% mycontrast %&% "$"
myfiles <- grep(myline, myfiles, value = TRUE)
myfiles <- paste0(inputFolder, "/", myfiles)

mygenes_group <- NULL
for (i in seq_len(length(myfiles))) {
  mygenes <- readRDS(myfiles[i])
  mygenes <- mygenes[order(mygenes$p_val_adj), ]
  if(!is.null(mypval)){
    mygenes <- mygenes[which(mygenes$p_val_adj < mypval),]
  }
  if(nrow(mygenes) >= mytop){
    mygenes_group <- c(mygenes_group, rownames(mygenes)[1:mytop])
  }else{
    mygenes_group <- c(mygenes_group, rownames(mygenes))
  }
}
mygenes_group <- unique(mygenes_group)

tmp <- paste(seurat_subset@meta.data[, which(colnames(seurat_subset@meta.data) == splitFeature)], seurat_subset@meta.data[, which(colnames(seurat_subset@meta.data) == contrastFeature)])
seurat_subset$myanno <- tmp


if (isLogNorm) {
  myexpression <- AverageExpression(seurat_subset, assays = myassay, features = mygenes_group, group.by = "myanno", slot = "data")
} else {
  myexpression <- AverageExpression(seurat_subset, assays = myassay, features = mygenes_group, group.by = "myanno", slot = "counts")
}
myexpression <- myexpression[[1]]

library(RColorBrewer)
hmcol <- rev(colorRampPalette(brewer.pal(8, "RdBu"))(127))
this_annot <- data.frame(rep(NA, ncol(myexpression)))
colnames(this_annot) <- splitFeature
rownames(this_annot) <- colnames(myexpression)
this_annot[, 1] <- unlist(lapply(strsplit(rownames(this_annot), " "), function(x) x[1]))
this_annot[, 2] <- unlist(lapply(strsplit(rownames(this_annot), " "), function(x) x[2]))
colnames(this_annot)[2] <- contrastFeature

phm <- pheatmap(t(myexpression),
  color = hmcol, scale = "column", fontsize = 10, fontsize_col = 7, fontsize_row = 7,
  show_colnames = TRUE, show_rownames = FALSE, annotation_names_col = TRUE,
  annotation_row = this_annot, # annotation_colors = annot_colors,
  main = "Top differentially expressed markers average normalized expression"
)
outname <- outdir %&% "DE_heatmap_" %&% splitFeature %&% "_" %&% contrastFeature %&% "_" %&% myassay %&% "_top" %&% mytop %&% "genes-" %&% mycontrast %&% ".png"
ggsave(outname, phm$gtable,
  width = 45, height = 20, units = "cm"
)
