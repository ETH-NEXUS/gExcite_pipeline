#################################################
## File name: integrateSamples_analyseCohort_ADT_GEX.R
## Author: Linda Grob
## Based on: integrateSamples_analyseCohort.R
## Date created: September 2020
## R Version: 3.5.1
##################################################

## GENERAL:
## This script is used for single cell integration analyses having ADT and GEX data using the Seurat R package.


## INPUT:
## - script takes list of RDS files (that contain SCE objects containing results of an
##   scRNA-seq experiment) and integrates those samples into one integrated Seurat object
##   (using the R package Seurat).
## - script takes second list of RDS files that contain results of GSVA analyses of
##   the same samples. IMPORTANT: Both lists have to be in the same order.
## - script can optionally work with an additional table of meta data about the samples.
## --> with package optparse for parsing command line arguments no parameter is strictly needed. If there is no meta data table to add simply do not use --add_meta
## For more details regarding input see help phrases of command line arguments.

##   META DATA FILE:
## Example:
##   sampleID     | cancerStage | pctRibosomalReads | reference
##   ----------------------------------------------------------
##   sample_1.RDS | stage III   | 0.3               | reference
##   sample_2.RDS | stage IV    | 0.1               | other
## If a plot should be generated from a meta data column the respective code block can be added in the Sample Visualisation section of the script.
## The Seurat function DimPlot() is used for categorical data and the function FeaturePlot() for continuous data.
## Note that DimPlot() does not plot a legend title using default settings (see additional line in existing code).

## OUTPUT:
## - integrated Seurat object in an RDS file.
## - UMAP plots showing GSVA results
## - UMAP plots characterising the cohort of samples. The plots include all integrated cells in one UMAP
##   and visualise data like cell types or g2m_score, n_umi, fractionMT.
## - UMAP plots showing gene expression
## - mean expression table per seurat cluster
## - UMAP showing combined adt & gex expression

## REMARK sampleID:
## The default is that the whole file name is entered as the sampleID (automatically retrieved from the input files and entered by the user in the meta data table).
## This is to keep this script very general and functioning without initial changes. But, especially for the plot
## showing the UMAP with all cells coloured by sampleID it is not optimal. If a shorter sampleID is desired you find example code with
## regular expressions to get a substring of the file name and have it serve as the 'sampleID'. This regular expression has to be changed and activated in two places
## indicated with 1) and 2).

## TODO: Infos about threshold and lookup
## TODO: Remove GSVA part
suppressPackageStartupMessages({
  library(Seurat)
  library(plyr)
  library(dplyr)
  library(SingleCellExperiment)
  library(optparse)
  library(reticulate)
  library(ggplot2)
  library(reshape2)
  library(cowplot)
  library(Matrix)
  library(patchwork)
  library(tidyverse)
})

# parse command line arguments
 option_list <- list(
   make_option("--cohort_list", type = "character", help = "List with paths to the input RDS files that contain SCE objects with GEX results. INPORTANT: GEX and GSVA input files have to be in the same order."),
 #  make_option("--gsva_list", type = "character", help = "List with paths to the input RDS files that contain SCE objects with GSVA results. IMPORTANT: GEX and GSVA input files have to be in the same order."),
   make_option("--colour_config", type = "character", help = "Path to config file table that sorts colours to cell types"),
 #  make_option("--selectedGenes", type = "character", help = "Tab delimited text file with the genes of interest in the first column. No header. The genes will be plotted and given out in chunks, alphabetical order."),
   make_option("--add_meta", type = "character", help = "Tab delimited text file with additional meta data. The first column must have header 'sampleID'. For a default run the full file name of the respective SCE analysis object is required, refer to script documentation on different sampleID formatting. Column headers describe categories. One line per input sample is required. Meta data file itself is optional: without table simply do not use --add_meta"),
#   make_option("--lookup", type = "character", help = "Full path to a lookup table matching protein to antibody names."),
#   make_option("--thresholds", type = "character", help = "Full path to a table listing adt thresholds for every sample included in the annalysis."),
   make_option("--outdir", type = "character", help = "Full path to output directory."),
   make_option("--outName", type = "character", help = "Prefix name of output files.")
 )
 opt_parser <- OptionParser(option_list = option_list)
 opt <- parse_args(opt_parser)


# convenience function for string concatenation
"%&%" <- function(a, b) paste(a, b, sep = "")
empty_string <- ""


########################
###   Read in data   ###
########################
# read in color config
config <- read.csv(opt$colour_config, sep = "\t", stringsAsFactors = FALSE)
print(config)

# read in data from all samples of the cohort
cohort_files <- read.csv(opt$cohort_list, header = FALSE, stringsAsFactors = FALSE)
all_sce_files <- cohort_files$V1
# make a list of all SCE objects
all_sce_objects <- lapply(seq_along(all_sce_files), function(x) {
  print(paste0("Read in sample ", x))
  print(all_sce_files[x])
  my_sce <- readRDS(all_sce_files[x])
  print(str(my_sce))
  # 1) example regex expression to shorten the sample name
  # example sample name is the basename of the input file up to (but excluding) the second "."
  # e.g. IDAMEME-T_scR_Ar1v7.6.genes_cells_filtered.corrected.atypical_removed.RDS --> IDAMEME-T_scR_Ar1v7.6
  my_sce@metadata$sample_name <- regmatches(basename(all_sce_files[x]), regexpr("[^\\.]+[^\\.]+", basename(all_sce_files[x])))
  #my_sce@metadata$sample_name <- basename(all_sce_files[x])
  # default sample name is the basename of the input file
  #my_sce@metadata$sample_name <- basename(all_sce_files[x])
  print(my_sce@metadata$sample_name)
  my_sce
})

print(empty_string)
# have number of samples
number_samples <- length(all_sce_files)
# sample name is the basename of the input file up to (but excluding) the second "."
all_sampleIDs <- sapply(all_sce_files, function(x) regmatches(basename(x), regexpr("[^\\.]+[^\\.]+", basename(x))))
# sample name is the basename of the input file
#all_sampleIDs <- sapply(all_sce_files, basename, simplify = TRUE)


# optionally, read in an additional meta data table
print(empty_string)
 if (is.character(opt$add_meta)) {
   addMeta_table <- read.csv(opt$add_meta, sep = "\t", stringsAsFactors = FALSE)
   ## Use sampleID up to (but excluding) the second "." as shorter sampleID
   addMeta_table$sampleID <- regmatches(addMeta_table$sampleID, regexpr("[^\\.]+[^\\.]+", addMeta_table$sampleID))
   print("Reading in meta data table: check if input samples and meta data sampleIDs are identical")
   stopifnot(sort(addMeta_table$sampleID) == sort(all_sampleIDs))
 } else if (is.null(opt$add_meta)) {
   print("No input file with additional meta data given.")
   addMeta_table <- NULL
 } else {
   print("Something went wrong reading in the meta data table.")
 }

###################################
###   Generate Seurat objects   ###
###################################

# Using Seurat SCTransform for sample integration
# Following the Vignette https://satijalab.org/seurat/v3.1/integration.html#SCTransform

# generate a list of Seurat objects from the list of SCE objects
# Fill the Seurat objects with the raw counts, gsva results, and metadata
print(empty_string)
all_samples <- lapply(seq_along(all_sce_objects), function(x) {
  print(paste0("Build seurat object of ", x, ". SCE object:"))
  # make sure that no cell and no gene (feature) is lost to default thresholds
  seurat_obj <- CreateSeuratObject(counts = assay(x = all_sce_objects[[x]], "counts"), min.cells = 0, min.features = 0)
  # create adt assay
  # get matrix containing only ADT counts
  # somehow R turns hyphens in column headers into dots (when using as.data.frame())
  # Therefore all hyphens are replaced by underscores.
  names_without_hyphen <- gsub("-", "_", names(colData(all_sce_objects[[x]])))
  colData_df <- as.data.frame(colData(all_sce_objects[[x]]))
  names(colData_df) <- names_without_hyphen
  mask_ADT_cols_temp <- names(colData_df)
  # Rename the col ADT_barcodes, to avoid including it latere on.
  mask_ADT_cols_temp <- replace(mask_ADT_cols_temp, mask_ADT_cols_temp=="ADT_barcodes", "Barcodes_ADT")
  print(mask_ADT_cols_temp)
  mask_ADT_cols <- grepl("ADT_",  mask_ADT_cols_temp)
  print(mask_ADT_cols)
  print(names(colData_df)[mask_ADT_cols])
  mask_ADT_cols_including_HT <- names(colData_df)[mask_ADT_cols]
  
  mask_ADT_cols_noHT <- !grepl("Hashtag",  mask_ADT_cols_including_HT)
  print("mask_ADT_cols_noHT")
  print(mask_ADT_cols_noHT)
  print(mask_ADT_cols_including_HT[mask_ADT_cols_noHT])
  adt_counts <- subset(colData_df, select = mask_ADT_cols_including_HT[mask_ADT_cols_noHT])
  adt_features <- colnames(adt_counts)[-1]
  seurat_obj[["adt"]] <- CreateAssayObject(counts =t(as.matrix(adt_counts)))
  # include all the metadata from the SCE objects
  seurat_obj[["phenograph_clusters"]] <- colData(all_sce_objects[[x]])$phenograph_clusters
  seurat_obj[["celltype_major"]] <- colData(all_sce_objects[[x]])$celltype_major
  seurat_obj[["celltype_final"]] <- colData(all_sce_objects[[x]])$celltype_final
  seurat_obj[["fractionMT"]] <- colData(all_sce_objects[[x]])$fractionMT
  seurat_obj[["n_umi"]] <- colData(all_sce_objects[[x]])$n_umi
  seurat_obj[["n_gene"]] <- colData(all_sce_objects[[x]])$n_gene
  seurat_obj[["log_umi"]] <- colData(all_sce_objects[[x]])$log_umi
  seurat_obj[["g2m_score"]] <- colData(all_sce_objects[[x]])$g2m_score
  seurat_obj[["s_score"]] <- colData(all_sce_objects[[x]])$s_score
  seurat_obj[["cycle_phase"]] <- colData(all_sce_objects[[x]])$cycle_phase
  seurat_obj[["celltype_major_full_ct_name"]] <- colData(all_sce_objects[[x]])$celltype_major_full_ct_name
  seurat_obj[["celltype_final_full_ct_name"]] <- colData(all_sce_objects[[x]])$celltype_final_full_ct_name
  seurat_obj[["sample_name"]] <- all_sce_objects[[x]]@metadata$sample_name
  seurat_obj[["phenograph_cluster"]] <- colData(all_sce_objects[[x]])$phenograph_clusters
  # add additional meta data from input table
  if (is.null(addMeta_table)) {
    print("No additional meta data is read into the seurat objects.")
  } else if (is.data.frame(addMeta_table)) {
    subset_meta <- subset(addMeta_table, sampleID == all_sce_objects[[x]]@metadata$sample_name)
    print("Check if subset of the meta data table only contains one sample.")
    stopifnot(dim(subset_meta)[1] == 1)
    # for each column of meta data add a column to the seurat object's meta data
    for (column in names(subset_meta)) {
      print(paste0("Adding column ", column, " to seurat object meta data."))
      seurat_obj[[column]] <- subset_meta[[column]]
    }
  } else {
    print("Something went wrong with the additional meta data table.")
  }
  # give out seurat object
  seurat_obj
})

##################################
###   Integration of samples   ###
##################################

# Perform SCTransform on each of the new Seuart objects, separately:
print(empty_string)
print("Perform SCTransform.")
print(Sys.time())
for (i in seq_len(length(all_samples))) {
  print("sctransform rna")
  all_samples[[i]] <- SCTransform(all_samples[[i]], assay = "RNA", verbose = FALSE, return.only.var.genes=FALSE)
}


# Perform ADT Normalization on each of the ADT assay, seperately:
for (i in seq_len(length(all_samples))) {
  #DefaultAssay(all_samples[[i]]) <- "adt"
  #RowsNA <- names(which(rowSums(is.na(all_samples[[i]]@assays$adt@counts))>0))
  #print(RowsNA)
  print("normalize adt")
  print(all_samples[[i]])
  all_samples[[i]] <- NormalizeData(all_samples[[i]], assay = "adt", normalization.method = "CLR", margin = 2)
  print("scale")
  all_samples[[i]] <- ScaleData(all_samples[[i]], assay = "adt")
  print("VariableFeatures")
  all_samples[[i]] <- FindVariableFeatures(all_samples[[i]], assay = "adt", selection.method = "vst") 
}


# Select integration features and prepare integration
print(empty_string)
print("Prepare integration GEX.")
print(Sys.time())

###################################################################################################################
###     Data integration with reference     ###
### data can be integrated using a reference data set. This reduces substantially the resources that are needed
### and can be the only option if large data sets are integrated, or if among the samples there are two or more
### samples that do not have any overlap in cell types.
### In this case the additional input meta data table must contain a column "reference" where the reference samples
### have the string "reference" and the others e.g. "other".
### Comment in the following lines and use the parameter "reference" in the function FindIntegrationAnchors
###################################################################################################################
# define reference data set
#reference_datasets <- which(unlist(lapply(all_samples, function(x) sum(x[[]]$reference == "reference") > 0)))
#print("Reference datasets:")
#print(reference_datasets)
#print(all_samples[reference_datasets])

# The number of features might have to be adapted per project, e.g. reduced if many samples are integrated or if no reference is used
integration_features <- SelectIntegrationFeatures(object.list = all_samples, nfeatures = 3000)
#integration_features <- SelectIntegrationFeatures(object.list = all_samples, nfeatures = 3000, assay = c("adt_gex","adt_gex","adt_gex"))
print("The following command yields in a warning if not all features are detected samples ")
all_samples <- PrepSCTIntegration(object.list = all_samples, anchor.features = integration_features)
# Find integration anchors
print(empty_string)
print("Find integration anchors GEX.")
print(Sys.time())

integration_anchors <- FindIntegrationAnchors(object.list = all_samples, normalization.method = "SCT",
                                  anchor.features = integration_features)
				#, reference = reference_datasets)

#integration_anchors <- FindIntegrationAnchors(object.list = all_samples, anchor.features = integration_features)
# Integrate samples
print(empty_string)
print("Integrate samples GEX.")
print(Sys.time())
seurat_integrated <- IntegrateData(anchorset = integration_anchors, normalization.method = "SCT", new.assay.name = "integratedGEX", integration_features)

# integrate ADT
for (i in seq_len(length(all_samples))) {
DefaultAssay(object = all_samples[[i]]) <- "adt"
}

print("integrate ADT")
print("Find integration anchors ADT")
integration_anchors <- FindIntegrationAnchors(object.list = all_samples)#,
                 #                 anchor.features = adt_features)
print("Integrate samples ADT")
seurat_integrated_temp <- IntegrateData(anchorset = integration_anchors, new.assay.name = "integratedADT")
print(str(seurat_integrated_temp))

integratedADT_data <- GetAssayData(seurat_integrated_temp, slot = "data", assay = "integratedADT")
dim(integratedADT_data)

# The following couple of lines have become necessary due to a Seurat Bug. 
# Could probably be written clearer once bug is resolved.
mat <- as.matrix(integratedADT_data)
mat <- mat[,colnames(seurat_integrated)]
features <- rownames(mat)
seurat_integrated[["integratedADT"]] <- CreateAssayObject(data= mat)
CreateAssayObject(data = mat)
print(str(seurat_integrated))

### Run normalisation (scale), PCA and UMAP, find neighbours and clusters on ADT.
print(empty_string)
print("Run normalisation (scale), PCA and UMAP, find neighbours and clusters on ADT.")
print(Sys.time())
DefaultAssay(object = seurat_integrated) <- "integratedADT"
seurat_integrated <- ScaleData(seurat_integrated, verbose = FALSE)
seurat_integrated <- RunPCA(object = seurat_integrated, verbose = TRUE, reduction.name = 'apca', features = features)

seurat_integrated <- RunUMAP(object = seurat_integrated, reduction = "apca", dims = 1:40)
seurat_integrated <- FindNeighbors(object = seurat_integrated, reduction = "apca")
seurat_integrated <- FindClusters(seurat_integrated, dims.use = 1:40)
seurat_integrated[["barcodes"]] <- rownames(seurat_integrated[[]])

print(str(seurat_integrated))
#seurat_integrated[["integratedADTscaled"]] <- CreateAssayObject(counts = GetAssayData(seurat_integrated_adt, slot="data"))

### if you need the integrated counts for more genes than the usual 1000-3000 that are the "integration_features"
# get feature/genes that will be integrated
# intersect of genes in cohort will be integrated
#to_integrate <- Reduce(intersect, lapply(integration.anchors@object.list, rownames))
# write into file all genes that will be integrated
#filename_integ <- paste0(opt$outdir, opt$outName, ".integrated_genes_intersect.tsv")
#write.table(to_integrate, filename_integ, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Union of all genes in the cohort
#to_integrate_all <- Reduce(union, lapply(integration.anchors@object.list, rownames))
#filename_integ_all <- paste0(opt$outdir, opt$outName, ".all_genes_in_cohort.tsv")
#write.table(to_integrate_all, filename_integ_all, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# run normalisation (scale), PCA and UMAP, find neighbours and clusters

print(empty_string)
print("Run normalisation (scale), PCA and UMAP, find neighbours and clusters on GEX.")
print(Sys.time())
DefaultAssay(object = seurat_integrated) <- "integratedGEX"
seurat_integrated <- ScaleData(seurat_integrated, verbose = FALSE)
seurat_integrated <- RunPCA(object = seurat_integrated,  verbose = TRUE, reduction.name = 'pca')

seurat_integrated <- RunUMAP(object = seurat_integrated, reduction = "pca", dims = 1:40)
seurat_integrated <- FindNeighbors(object = seurat_integrated)
seurat_integrated <- FindClusters(seurat_integrated, dims.use = 1:40)
seurat_integrated[["barcodes"]] <- rownames(seurat_integrated[[]])


# Find multimodal neighbors
seurat_integrated <- FindMultiModalNeighbors(
  seurat_integrated, reduction.list = list("pca", "apca"), 
  dims.list = list(1:50, 1:50), modality.weight.name = "RNA.weight"
)


seurat_integrated <- RunUMAP(seurat_integrated, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
seurat_integrated <- FindClusters(seurat_integrated, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)


seurat_integrated <- RunUMAP(seurat_integrated, reduction = 'pca', dims = 1:50, assay = 'RNA', 
              reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
seurat_integrated <- RunUMAP(seurat_integrated, reduction = 'apca', dims = 1:50, assay = 'ADT', 
              reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')



###############################
###   Sample Visualisation  ###
###############################
DefaultAssay(object = seurat_integrated) <- "integratedGEX"
# read in colour_config and sort colours to cell types
config$cell_type <- gsub(pattern = "([^_]+)_.*", "\\1", config$cell_type)
print("Check if the number of cell types and the number of colours for the cell types match.")
stopifnot(length(config$colour) == length(config$cell_type))
stopifnot(length(unique(config$colour)) == length(unique(config$mapping)))
config$class <- NULL
config$indication <- NULL

# have character vector of colour names from the colour config

config <- rbind(config, c("uncertain", "grey50", "uncertain"))
config <- rbind(config, c("unknown", "black", "unknown"))
# give names to colours from other column of colour config
ct.color <- config$colour
names(ct.color) <- config$mapping
# make sure 'celltype_major' and 'celltype_final' metadata columns contain factors and are mapped to the correct celltype
seurat_integrated[["celltype_major"]] <- as.factor(config$mapping[match(as.factor(seurat_integrated[[]]$celltype_major), config$cell_type)])
seurat_integrated[["celltype_final"]] <- as.factor(config$mapping[match(as.factor(seurat_integrated[[]]$celltype_final), config$cell_type)])

# make sure "original objects" contain the correct mapping for the celltypes
for (i in seq_len(length(all_samples))) {
  #print(which(is.na(all_samples[[i]]$celltype_final )))
  all_samples[[i]]$celltype_major <- as.factor(config$mapping[match(as.factor(all_samples[[i]]$celltype_major), config$cell_type)])
  all_samples[[i]]$celltype_final <- as.factor(config$mapping[match(as.factor(all_samples[[i]]$celltype_final), config$cell_type)])
}

# Write integrated Seurat object into RDS file
 filename_out <- opt$outdir %&% opt$outName %&% ".integrated_seurat_object.RDS"
 print(paste0("File name of output RDS file with integrated Seurat object: ", filename_out))
 saveRDS(seurat_integrated, filename_out)
