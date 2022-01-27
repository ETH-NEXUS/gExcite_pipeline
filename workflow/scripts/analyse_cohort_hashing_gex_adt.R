###################################################
## File name: analyse_cohort_hashing_gex_adt.R
## Author: Anne Bertolini, Linda Grob
## Date created: August 2020
## Extended: October 2020
## R Version: 3.5.1
###################################################

## GENERAL:
## This script combines GEX and ADT results of different samples into one object.
## Those results can then be visualised together.
## No integration or normalisation step is included.


library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(scMerge)
library(uwot)
library(rhdf5)
library(stringr)
library(patchwork)
library(optparse)
library(limma)
library(cowplot)
library(reshape2)
library(ggridges)
library(igraph)
library(BREMSC)

# parse command line arguments
option_list <- list(
  make_option("--hashing_results", type = "character", help = "Table with results from demultiplexing of hashed samples."),
  make_option("--metadata_table", type = "character", help = "Table with columns RDS_files (paths to input RDS files that contain SCE objects with GEX and ADT results), h5_files (paths to h5 files that contain GEX most variable genes), hashtag (hashtag corresponding to the RDS file, . if not available) and sample_name."),
  make_option("--colour_config", type = "character", help = "Colour config for visualising the cell types. Sometimes union of different cell type colour config files is required."),
  make_option("--lookup", type = "character", help = "Path to a lookup table for gene & protein names."),
  make_option("--threshold", type = "character", help = "File containing individual log thresholds for each antibody."),
  make_option("--sampleName", type = "character", help = "Sample name that will be prefix of all plot names."),
  make_option("--outdir", type = "character", help = "Path to the outdir, in which files will be written."),
  make_option("--variableGenesSample", type = "integer", help = "Number of HVGs that is used to reduce the pearson residuals matrix from the gex data.", default = 500),
  make_option("--threads", type = "integer", help = "Number of available threads to run the script.", default = 1)
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# convenience function for string concatenation
"%&%" <- function(a, b) paste(a, b, sep = "")
empty_string <- ""


###################################
### General plotting parameters ###
###################################

# aspect ratio of UMAP plots
UMAP_aspect_ratio <- "1"
# chunk size of gene sets in GSVA plots
chunksize_genesets <- 12
# number of columns in GSVA plots
columns_genesets <- 4


###################################
###   Read in and format data   ###
###################################

# read in metadata table
cat("\nRead in metadata table:\n")
metadata <- read.table(opt$metadata_table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
print(metadata)

# read in hashtags
hash <- metadata$hashtag
no_tags <- grepl("\\.", hash)
print(no_tags)
true_tags <- hash[!no_tags]
print(true_tags)
exclude_hashplots <- FALSE
# Exclude the plots regarding the hashing analysis if not all samples where hashed together.
if (length(true_tags) < length(hash)) {
  exclude_hashplots <- TRUE
}
print(exclude_hashplots)
# read in data from all samples of the cohort
cat("\nRead in one SCE object per sample:\n")
all_sce_files <- metadata$RDS_files
cell_desc <- character()

# get union of samples' most variable genes in the GEX pipeline
# only those will be used for calculation of UMAP
var_genes_files <- metadata$h5_files

var_genes_list <- lapply(seq_along(var_genes_files), function(x) {
  cat("\n\nRead in variable genes sample ", x, "\n")
  cat(var_genes_files[x], "\n")
  var_genes_temp <- rhdf5::h5read(var_genes_files[x], "gene_attrs/gene_names")
  var_genes_temp
})
var_genes <- unique(unlist(var_genes_list))

# make a list of all SCE objects
all_sce_objects <- lapply(seq_along(all_sce_files), function(x) {
  cat("\n\nRead in sample ", x, "\n")
  cat(all_sce_files[x], "\n")
  my_sce <- readRDS(all_sce_files[x])
  # default sample name is the basename of the input file
  my_sce@metadata$sample_name <- metadata$sample_name[x]
  cat(my_sce@metadata$sample_name, "\n")
  colData(my_sce)$sampleID <- metadata$sample_name[x]
  colData(my_sce)$file_name <- basename(all_sce_files[x])
  my_sce
})

# make a list of all gsva SCE objects
all_gsva_files <- metadata$gsva
cat("Check if number of single cell pipeline results and GSVA results are the same.\n")
stopifnot(length(all_gsva_files) == length(all_sce_files))
# make a list of all GSVA results that are stored in SCE objects as well
all_gsva_results <- lapply(seq_along(all_gsva_files), function(x) {
  cat("Read in GSVA results ", x, "\n")
  cat(all_gsva_files[x], "\n")
  my_gsva <- readRDS(all_gsva_files[x])
  # gsva sample name is the basename of the gsva input file. Name is only used for printing out per sample pair, as sanity check.
  # Matching is solely based on the order of the files.
  my_gsva@metadata$sample_name <- metadata$sample_name[x]
  my_gsva@metadata$file_name <- basename(all_gsva_files[x])
  print(my_gsva@metadata$sample_name)
  stopifnot(colnames(my_gsva) == colData(my_gsva)$barcodes)
  stopifnot(rownames(colData(my_gsva)) == colData(my_gsva)$barcodes)
  my_gsva
})

# check if objects in all_sce_objects and all_gsva_results have same order
# and contain the same cells
# add gsva results to the SCE object with GEX data (to colData)
cat("Checking if sce_object and gsva_result contain the same cells\n\n")
for (i in seq_along(all_sce_files)) {
  cat("List entry number", i, "), see names of the sample pair (GEX and GSVA):\n")
  cat(all_sce_objects[[i]]@metadata$sample_name, "\n")
  cat(all_gsva_results[[i]]@metadata$sample_name, "\n\n")
  cat("Check if the cells in gex SCE object are the same as in GSVA SCE object.\n")
  barcodes_sce <- colnames(all_sce_objects[[i]])
  all_gsva_results[[i]] <- all_gsva_results[[i]][, barcodes_sce]
  stopifnot(colnames(all_sce_objects[[i]]) == colnames(all_gsva_results[[i]]))
  # data frame with gsva results, cells are rows
  gsva_df <- as.data.frame(t(assay(all_gsva_results[[i]], "gsva")))
  # check if cell barcodes are the same
  stopifnot(colnames(all_sce_objects[[i]]) == rownames(gsva_df))
  # Add"GSVA_" to GSVA headers to retrieve them later
  names(gsva_df) <- paste0("GSVA_", names(gsva_df))
  # add all gsva results to the colData table
  for (gsva_column in names(gsva_df)) {
    colData(all_sce_objects[[i]])[gsva_column] <- gsva_df[[gsva_column]]
  }
}

# make a list of all highest variable genes based on the number in opt$variableGenesSample
print("Filtering the highest variable genes down using the opt$variableGenesSample")
hvg_list <- lapply(seq_along(all_sce_objects), function(x) {
  subset_vari_rowData <- rowData(all_sce_objects[[x]])[rownames(rowData(all_sce_objects[[x]])) %in% var_genes, ]
  subset_vari_ordered <- subset_vari_rowData[order(subset_vari_rowData$residual_variance, decreasing = TRUE), ]
  subset_vari_ordered$gene_names[1:opt$variableGenesSample]
})
var_genes <- unique(unlist(hvg_list))
print("Number of highest variable genes finally used:")
print(length(var_genes))

# read in hashing results table
if (exclude_hashplots == FALSE) {
  cat("\nRead in table with hashing results:\n")
  cat("\n", opt$hashing_results, "\n")
  hashing_results <- read.csv(opt$hashing_results, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}

# read in lookup table
cat("\nRead in lookup table that connects ADTs to genes of interest:\n")
cat("\n", opt$lookup, "\n")
lookup <- read.csv(opt$lookup, header = TRUE, sep = "\t")
lookup$Antibody <- gsub("-", "_", lookup$Antibody)
cat("\nThe following lookup table is used to map gene names to antibodies and the other way round.\n")
print(lookup)

# read in threshold file
print("Read in threshold file.")
thresholds <- read.table(opt$threshold, header = TRUE, sep = "\t", row.names = 1)
names_without_hyphen <- gsub("-", "_", rownames(thresholds))
rownames(thresholds) <- names_without_hyphen
print("Checking if for every sampleName a threshold column is given.")
print(metadata$sample_name)
stopifnot(all(metadata$sample_name %in% colnames(thresholds)))
print("The following thresholds have been read for this sample:")
print(thresholds)

print(
  "Script assumes input thresholds to be logged. Raw Thresholds are then calculated to be the following:"
)
thresholds_exp <- round(exp(thresholds))
print(thresholds_exp)
# read in colour config
cat("\nRead in Color-Config\n")
cat("\n", opt$colour_config, "\n")
config <- read.csv(opt$colour_config, sep = "\t", stringsAsFactors = FALSE)
config$cell_type <- gsub(pattern = "([^_]+)_.*", "\\1", config$cell_type)
stopifnot(length(config$colour) == length(config$cell_type))
stopifnot(config$colour == unique(config$colour))
ct.color <- c(config$colour, "grey50", "black")
names(ct.color) <- c(config$cell_type, "uncertain", "unknown")
print(ct.color)

# get list of all cell description columns (colData) in all samples to check if
# GEX cell description data is identical between samples
cell_desc_list <- lapply(seq_along(all_sce_objects), function(x) {
  names(colData(all_sce_objects[[x]]))
})
cell_desc <- unlist(cell_desc_list)
# get unique list
unique_cell_desc <- unique(cell_desc)
# get frequencies of colData headers. Should be identical, and always number of samples.
freq_cell_desc <- as.data.frame(table(cell_desc))
unique(freq_cell_desc$Freq)
cat("\n\nCheck if the cell description columns in all samples are identical.\n")
stopifnot(unique(freq_cell_desc$Freq) == length(cell_desc_list))


###   combine all SCE objects   ###
combined <- scMerge::sce_cbind(all_sce_objects,
  exprs = c("counts", "normcounts", "pearson_resid"),
  method = "union",
  cut_off_batch = 0.0,
  cut_off_overall = 0.0,
  colData_names = unique_cell_desc
)

# find duplicated barcodes
barcodes_no_prefix <- gsub(pattern = "[^-]+-(.*)", "\\1", colnames(combined))
indices_duplicated <- which(duplicated(barcodes_no_prefix))
barcodes_duplicated <- barcodes_no_prefix[indices_duplicated]
indices_ALL_duplicated <- which(barcodes_no_prefix %in% barcodes_duplicated)
mask_ALL_duplicated <- barcodes_no_prefix %in% barcodes_duplicated
cat("\n\nDuplicated barcodes:\n")
print(table(mask_ALL_duplicated))
cat("\n")

# Remove all duplicated barcodes
combined <- combined[, !mask_ALL_duplicated]

# add hashing result table to SCE object
if (exclude_hashplots == FALSE) {
  cat("Add hashing results to SCE object colData")
  # get subset of hashing results containing only cells that are in SCE object
  hashing_results_subset <- hashing_results[hashing_results$barcodes %in% colData(combined)$barcodes, ]
  barcodes_intersect <- intersect(colData(combined)$barcodes, hashing_results_subset$barcodes)
  stopifnot(length(barcodes_intersect) == length(hashing_results_subset$barcodes))
  stopifnot(length(barcodes_intersect) == length(colData(combined)$barcodes))
  # make sure the hashing data is in the right order
  ordered_barcodes <- as.data.frame(colData(combined)$barcodes)
  names(ordered_barcodes) <- "barcodes"
  hashing_results_subset <- plyr::join(ordered_barcodes, hashing_results_subset, by = "barcodes")
  stopifnot(colData(combined)$barcodes == hashing_results_subset$barcodes)
  # add each column of the hashing data to the SCE object colData
  for (column in names(hashing_results_subset)) {
    cat("\nAdding column", column, "from hashing table to combined SCE object.\n")
    colData(combined)[column] <- hashing_results_subset[[column]]
  }
}

###################################
###        Ridgeplot ABs        ###
###################################

# Normalise ADT counts using thresholds or simply substract thresholds
# get matrix containing only ADT counts
# somehow R turns hyphens in column headers into dots (when using as.data.frame())
# Therefore all hyphens are replaced by underscores.
names_without_hyphen <- gsub("-", "_", names(colData(combined)))
colData_df <- as.data.frame(colData(combined))
names(colData_df) <- names_without_hyphen
# extract sampleID together with ADT columns
mask_ADT_cols <- grepl("sampleID|ADT_", names(colData_df))
adt_counts <- subset(colData_df, select = names(colData_df)[mask_ADT_cols])
adt_counts$ADT_barcodes <- NULL

# Currently here only values below threshold are set to NA and later on coloured gray.
# Nevertheless it makes probably sense to keep the function for the time being, since we might want to normalize later.
normalize_adt <- function(x, threshold) {
  x[x < threshold] <- NA
  #   return ((y - quantile(y, 0.01)) / (quantile(y, 0.99) - quantile(y, 0.01)))
  return(x)
}

adt_only_names <- names(adt_counts)[grepl("ADT_", names(adt_counts))]
for (sampleID in metadata$sample_name) {
  for (ab in adt_only_names) {
    ab_threshold <- thresholds[substring(ab, 5), sampleID]
    subset <- adt_counts[adt_counts$sampleID == sampleID, ab]
    normalized_subset <- normalize_adt(subset, ab_threshold)
    # Cap the values to the interval [0,1]
    #       normalized_subset[normalized_subset < 0] <- 0
    #      normalized_subset[normalized_subset > 1] <- 1
    adt_counts[adt_counts$sampleID == sampleID, ab] <- normalized_subset
  }
}

adt_counts$sampleID <- NULL
names(adt_counts) <- paste("Norm_", names(adt_counts), sep = "")

# add each column of the norm_ADT data to the SCE object colData
for (column in names(adt_counts)) {
  cat("\nAdding column", column, "from adt_counts table to combined SCE object.\n")
  colData(combined)[column] <- adt_counts[[column]]
}

# get matrix containing only ADT counts
# somehow R turns hyphens in column headers into dots (when using as.data.frame())
# Therefore all hyphens are replaced by underscores.
names_without_hyphen <- gsub("-", "_", names(colData(combined)))
colData_df <- as.data.frame(colData(combined))
names(colData_df) <- names_without_hyphen
mask_ADT_cols <- grepl("^ADT_", names(colData_df))
adt_counts <- subset(colData_df, select = names(colData_df)[mask_ADT_cols])
adt_counts$ADT_barcodes <- NULL
adt_counts_mat <- as.matrix(adt_counts)
# adt counts without hashtags
mask_Hashtag_cols <- grepl("_Hashtag", names(adt_counts), ignore.case = TRUE)
adt_counts_no_hashtags <- subset(adt_counts, select = names(adt_counts)[!mask_Hashtag_cols])

# have matrix with GEX pearson residuals of all variable genes
residuals_variableGenes <- assay(combined, "pearson_resid")[rownames(combined) %in% var_genes, ]
# rows must contain observations (cells)
residuals_variableGenes <- t(residuals_variableGenes)


###########################
###   Calculate UMAPS   ###
###########################

# get four UMAPS, based on adt data, gex data, both (adt + gex), and finally adt without the hashtag counts
set.seed(3792)
for (type in c("adt", "gex", "adt_gex_withHashtags", "adt_noHashtags", "adt_gex")) {
  if (type == "adt") {
    cat("\nGet ADT based UMAP embedding.\n")
    print(names(adt_counts))
    # calculate UMAP coordinates
    umap_adt <- uwot::umap(adt_counts, n_neighbors = 30, pca = 50, spread = 1, min_dist = 0.1)
    # reformat UMAP coordinates
    umap_adt_df <- as.data.frame(umap_adt)
    names(umap_adt_df) <- c("ADT_umap1", "ADT_umap2")
    umap_adt_df$barcodes <- rownames(adt_counts)
    # make sure all cell dimensions are the same
    stopifnot(combined$barcodes == umap_adt_df$barcodes)
    # add UMAP coordinates to SCE object
    reducedDim(combined, "umap_adt") <- umap_adt_df
    # add umap coordinates to the meta data table of cells
    colData_df <- plyr::join(colData_df, umap_adt_df, by = "barcodes")
  } else if (type == "gex") {
    cat("\nGet GEX based UMAP embedding.\n")
    # TODO abgleichen mit der UMAP Erstellung im preprocessing. 2000 Most highly variable genes?
    # use union of all highly variable genes to calculate UMAP
    umap_gex <- uwot::umap(residuals_variableGenes,
      n_neighbors = 30, pca = 50, spread = 1, min_dist = 0.1
    )
    # reformat UMAP coordinates
    umap_gex_df <- as.data.frame(umap_gex)
    names(umap_gex_df) <- c("GEX_umap1", "GEX_umap2")
    umap_gex_df$barcodes <- rownames(residuals_variableGenes)
    # make sure all cell dimensions are the same
    stopifnot(combined$barcodes == umap_gex_df$barcodes)
    # add UMAP coordinates to SCE object
    reducedDim(combined, "umap_gex") <- umap_gex_df
    # add umap coordinates to the meta data table of cells
    colData_df <- plyr::join(colData_df, umap_gex_df, by = "barcodes")
  } else if (type == "adt_gex_withHashtags") {
    cat("\nGet ADT and GEX based UMAP embedding.\n")
    # TODO why is +1 added to the ADT counts?
    # adt_counts_plus1 <- adt_counts + 1
    # only keep cells that are have ADT counts and GEX pearson residuals
    adt_counts_overlapping <- adt_counts[rownames(adt_counts) %in% rownames(residuals_variableGenes), ]
    gex_residuals_overlapping <- residuals_variableGenes[rownames(residuals_variableGenes) %in% rownames(adt_counts), ]
    stopifnot(rownames(adt_counts) == rownames(residuals_variableGenes))
    # Combine ADT counts and GEX pearson residuals to one matrix
    ADT_GEX_combined <- cbind(adt_counts_overlapping, gex_residuals_overlapping)
    # calculate UMAP coordinates
    umap_adt_gex <- uwot::umap(ADT_GEX_combined,
      scale = TRUE, n_neighbors = 30, pca = 50, spread = 1, min_dist = 0.1, ret_nn = T
    )
    # reformat UMAP coordinates
    umap_adt_gex_df <- as.data.frame(umap_adt_gex$embedding)
    names(umap_adt_gex_df) <- c("ADT_GEX_withHashtags_umap1", "ADT_GEX_withHashtags_umap2")
    umap_adt_gex_df$barcodes <- rownames(ADT_GEX_combined)
    # add UMAP coordinates to SCE object
    reducedDim(combined, "umap_adt_gex_withHashtags") <- umap_adt_gex_df
    # add umap coordinates to the meta data table of cells
    colData_df <- plyr::join(colData_df, umap_adt_gex_df, by = "barcodes")
  } else if (type == "adt_gex") {
    cat("\nGet ADT and GEX based UMAP embedding without hashtag counts.\n")
    print(names(adt_counts_no_hashtags))
    # only keep cells that have ADT counts and GEX pearson residuals
    adt_counts_overlapping <- adt_counts_no_hashtags[rownames(adt_counts_no_hashtags) %in% rownames(residuals_variableGenes), ]
    gex_residuals_overlapping <- residuals_variableGenes[rownames(residuals_variableGenes) %in% rownames(adt_counts_no_hashtags), ]
    stopifnot(rownames(adt_counts) == rownames(residuals_variableGenes))
    # Combine ADT counts and GEX pearson residuals to one matrix
    ADT_GEX_combined <- cbind(adt_counts_overlapping, gex_residuals_overlapping)
    # calculate UMAP coordinates
    umap_adt_gex_noHashtags <- uwot::umap(ADT_GEX_combined,
      scale = TRUE, n_neighbors = 30, pca = 50, spread = 1, min_dist = 0.1, ret_nn = T
    )
    # reformat UMAP coordinates
    umap_adt_gex_noHashtags_df <- as.data.frame(umap_adt_gex_noHashtags$embedding)
    names(umap_adt_gex_noHashtags_df) <- c("ADT_GEX_umap1", "ADT_GEX_umap2")
    umap_adt_gex_noHashtags_df$barcodes <- rownames(ADT_GEX_combined)
    # add UMAP coordinates to SCE object
    reducedDim(combined, "umap_adt_gex") <- umap_adt_gex_noHashtags_df
    # add umap coordinates to the meta data table of cells
    colData_df <- plyr::join(colData_df, umap_adt_gex_noHashtags_df, by = "barcodes")
    } else if (type == "adt_noHashtags") {
    cat("\nGet ADT based UMAP embedding without hashtag counts.\n")
    print(names(adt_counts_no_hashtags))
    # calculate UMAP coordinates
    umap_adt_noHashtags <- uwot::umap(adt_counts_no_hashtags, n_neighbors = 30, pca = 50, spread = 1, min_dist = 0.1)
    # reformat UMAP coordinates
    umap_adt_noHashtags_df <- as.data.frame(umap_adt_noHashtags)
    names(umap_adt_noHashtags_df) <- c("ADT_noHashtags_umap1", "ADT_noHashtags_umap2")
    umap_adt_noHashtags_df$barcodes <- rownames(adt_counts_no_hashtags)
    # make sure all cell dimensions are the same
    stopifnot(combined$barcodes == umap_adt_noHashtags_df$barcodes)
    # add UMAP coordinates to SCE object
    reducedDim(combined, "umap_adt_noHashtags") <- umap_adt_noHashtags_df
    # add umap coordinates to the meta data table of cells
    colData_df <- plyr::join(colData_df, umap_adt_noHashtags_df, by = "barcodes")
  }
}


##########################
###   Generate Plots   ###
##########################

# match colours to cell type levels
id.final.ct <- match(levels(as.factor(colData_df$celltype_final)), names(ct.color))
id.major.ct <- match(levels(as.factor(colData_df$celltype_major)), names(ct.color))
celltype_list <- list(id.final.ct, id.major.ct)
names(celltype_list) <- c("celltype_final", "celltype_major")

# plot cell related data to characterize the cohort
# plot it on all three UMAPs
for (type in c("ADT", "GEX", "ADT_GEX_withHashtags", "ADT_noHashtags","ADT_GEX")) {
  cat("\n\nGenerate plots that characterise cohort. Use", type, "UMAP coordinates\n\n")
  umap1 <- paste0(type, "_umap1")
  umap2 <- paste0(type, "_umap2")

  outdir <- opt$outdir
  if (type %in% c("ADT", "GEX", "ADT_noHashtags","ADT_GEX_withHashtags")) {
    outdir <- paste0(opt$outdir, type, "_based_embedding/")
    dir.create(outdir)
  }
  # plot sample IDs ("batch")
  print("# plot sample IDs ")
  # generate sample UMAP without legend for later use in adt-gex overview plots
  umap_sample <- ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = sampleID)) +
    geom_point(size = 1) +
    xlab(umap1) +
    ylab(umap2) +
    theme_light(base_size = 15) +
    theme(aspect.ratio = UMAP_aspect_ratio) +
    theme(legend.position = "none")
  # generate sample UMAP with legend to save it seperatly
  umap_sample_final <- ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = sampleID)) +
    geom_point(size = 1) +
    xlab(umap1) +
    ylab(umap2) +
    theme_light(base_size = 15) +
    guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
    theme(aspect.ratio = UMAP_aspect_ratio)
  legend_sample <- cowplot::get_legend(umap_sample_final)
  filename <- paste0(outdir, opt$sampleName, "_", type, "__samples.png")
  ggsave(umap_sample_final, filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

  # plot file names (sanity check)
  umap_filename <- ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = file_name)) +
    geom_point(size = 1) +
    xlab(umap1) +
    ylab(umap2) +
    theme_light(base_size = 15) +
    guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
    theme(aspect.ratio = UMAP_aspect_ratio)
  filename <- paste0(outdir, opt$sampleName, "_", type, "__fileNames.png")
  ggsave(umap_filename, filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

  # plot cell cycle phase
  print("# plot cell cycle phase")
  ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = cycle_phase)) +
    geom_point(size = 1) +
    xlab(umap1) +
    ylab(umap2) +
    theme_light(base_size = 15) +
    guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
    theme(aspect.ratio = UMAP_aspect_ratio)
  filename <- paste0(outdir, opt$sampleName, "_", type, "__cellCycle.png")
  ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

  # plot first cell types
  print("# plot first cell type")
  ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = celltype_major)) +
    scale_color_manual(
      name = "Cell type",
      values = ct.color[id.major.ct],
      drop = F
    ) +
    geom_point(size = 1) +
    xlab(umap1) +
    ylab(umap2) +
    theme_light(base_size = 15) +
    guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 20)) +
    theme(
      aspect.ratio = UMAP_aspect_ratio,
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 14),
      legend.key.width = unit(0.5, "line"),
      legend.key.height = unit(0.5, "line"),
      legend.spacing.y = unit(0.5, "line"),
      legend.spacing.x = unit(0.5, "line")
    )
  filename <- paste0(outdir, opt$sampleName, "_", type, "__celltypes_major.png")
  ggsave(filename = filename, width = 35, height = 20, dpi = 300, units = "cm")

  # plot second cell types
  print("# plot second cell type")
  ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = celltype_final)) +
    scale_color_manual(
      name = "Cell type",
      values = ct.color[id.final.ct],
      drop = F
    ) +
    geom_point(size = 1) +
    xlab(umap1) +
    ylab(umap2) +
    theme_light(base_size = 15) +
    guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 20)) +
    theme(
      aspect.ratio = UMAP_aspect_ratio,
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 14),
      legend.key.width = unit(0.5, "line"),
      legend.key.height = unit(0.5, "line"),
      legend.spacing.y = unit(0.5, "line"),
      legend.spacing.x = unit(0.5, "line")
    )
  filename <- paste0(outdir, opt$sampleName, "_", type, "__celltypes_final.png")
  ggsave(filename = filename, width = 35, height = 20, dpi = 300, units = "cm")

  # take natural log of each hashtag (count + 1) into new column
  print("# take natural log of each hashtag (count + 1) into new column")

  for (hashtag in paste("ADT_", true_tags, sep = "")) {
    colData_df[[paste0("ln_", hashtag)]] <- round(log(colData_df[[hashtag]] + 1))
  }
  # set several continuous values to be plotted
  print("# set several continuous values to be plotted")
  values_to_plot <- c("fractionMT", "log_umi", "n_gene", paste("ln_ADT_", true_tags, sep = ""))
  if (exclude_hashplots == FALSE) {
    values_to_plot <- c("HTO_margin_1", "HTO_margin_2", "fractionMT", "log_umi", "n_gene", paste("ln_ADT_", true_tags, sep = ""))
  }
  print(values_to_plot)
  # plot various continuous values
  for (continuousValue in values_to_plot) {
    cat("\nPlotting", continuousValue, "\n")
    ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = colData_df[[continuousValue]])) +
      geom_point(size = 1) +
      xlab(umap1) +
      ylab(umap2) +
      theme_light(base_size = 15) +
      theme(aspect.ratio = UMAP_aspect_ratio) +
      scale_colour_viridis_c(name = continuousValue)
    filename <- paste0(outdir, opt$sampleName, "_", type, "__", continuousValue, ".png")
    ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")
  }
}


#######################################
###            GSVA plots           ###
#######################################
cat("\n\nGSVA plots\n")
print(Sys.time())
mask_GSVA_cols <- grepl("GSVA_", names(colData_df))
gsva_counts <- subset(colData_df, select = names(colData_df)[mask_GSVA_cols])
gsva_counts$barcodes <- colData_df$barcodes
# remove the "GSVA_" tag
colnames(gsva_counts) <- gsub(pattern = "[^_]+_(.*)", "\\1", colnames(gsva_counts))
gsva_melted <- reshape2::melt(gsva_counts, id.vars = "barcodes")

# Trim outlier values and have one scale for all
gsva_melted$value_trimmed <- gsva_melted$value
limit_min <- quantile(gsva_melted$value, prob = 0.01)
gsva_melted$value_trimmed[gsva_melted$value <= limit_min] <- limit_min
print("Number of data points set to limit_min:")
print(table(gsva_melted$value_trimmed == limit_min))
limit_max <- quantile(gsva_melted$value, prob = 0.99)
gsva_melted$value_trimmed[gsva_melted$value >= limit_max] <- limit_max
print("Number of data points set to limit_max:")
print(table(gsva_melted$value_trimmed == limit_max))
print("summary(gsva_melted$value_trimmed) after trimmming the outliers:")
print(summary(gsva_melted$value_trimmed))
# get breaks
gsva_melted$value_trimmed <- round(gsva_melted$value_trimmed, digits = 2)
break_low <- min(gsva_melted$value_trimmed)
break_high <- max(gsva_melted$value_trimmed)
label_low <- paste("<=", break_low)
print(label_low)
label_high <- paste(">=", break_high)
print(label_high)
print("summary(gsva_melted$value_trimmed) after rounding the values to two digits:")
print(summary(gsva_melted$value_trimmed))

# get chunks of gene sets
all_genesets <- colnames(gsva_counts)
chunks_genesets <- split(all_genesets, ceiling(seq_along(all_genesets) / chunksize_genesets))
# get umap coordinates of the integrated cells
umap_coordinates <- subset(colData_df, select = c("ADT_GEX_umap1", "ADT_GEX_umap2"))
umap_coordinates$barcodes <- colData_df$barcodes

# set plotting theme for gsva plots
theme_set(theme_grey())
# generate plot per chunk of gene sets
for (chunk in seq_along(chunks_genesets)) {
  cat("\n\nChunk: ", chunk, "\n")
  cat("Plot chunk of gene sets:\n")
  print(chunks_genesets[[chunk]])
  sub_gsva_melted <- gsva_melted[gsva_melted$variable %in% chunks_genesets[[chunk]], ]
  gsva_data_plotting <- plyr::join(sub_gsva_melted, umap_coordinates, type = "inner")
  cat("Check if the cell barcodes in both melted gsva objects are identical.\n")
  stopifnot(gsva_data_plotting$barcodes == sub_gsva_melted$barcodes)
  # plot
  p_genesets <- ggplot(gsva_data_plotting, aes(x = ADT_GEX_umap1, y = ADT_GEX_umap2)) +
    geom_point(aes(color = value_trimmed), size = 1) +
    theme(
      panel.background = element_rect(fill = "grey80"),
      strip.text.x = element_text(size = 10.5),
      aspect.ratio = UMAP_aspect_ratio
    ) +
    scale_color_distiller(
      name = "", palette = "RdBu",
      breaks = c(break_low, 0, break_high),
      labels = c(label_low, "0.00", label_high)
    ) +
    facet_wrap(~variable, ncol = columns_genesets)
  p_genesets
  filename_gsva_chunk <- paste0(opt$outdir, opt$sampleName, ".integrated_gsva_umap__", chunk, ".png")
  ggsave(
    filename = filename_gsva_chunk,
    width = 40, height = 7.5 * ceiling(length(chunks_genesets[[chunk]]) / columns_genesets),
    dpi = 300, units = "cm"
  )
}


#######################################
###   Antibody / Expression plots   ###
#######################################

# specify one UMAP version for the expression plots
type_exp <- "ADT_GEX"
cat("\n\nGenerate UMAP expression plots. Use", type_exp, "UMAP coordinates\n\n")
umap1_exp <- paste0(type_exp, "_umap1")
umap2_exp <- paste0(type_exp, "_umap2")


print("Start clustering:")
# Coulours used later on for clustering
cols33 <- c(
  "red2", "green4", "blue2", "cyan2", "yellow1", "purple", "brown",
  "chocolate1", "chartreuse2", "darkgoldenrod3", "steelblue1", "slateblue3", "olivedrab4", "gold2",
  "violetred3", "darkcyan", "orchid3", "darksalmon", "darkslategrey", "khaki", "indianred2", "magenta", "slategray2",
  "olivedrab1", "mediumaquamarine", "hotpink", "yellow3",
  "bisque4", "darkseagreen1", "dodgerblue3",
  "deeppink4", "sienna4", "mediumorchid4"
)
print("Defined colours for clustering; only 33 available, if more clusters are found, this might crash.")
### Cluster
print("Start clustering on ADT and GEX data.")
# In order to have an idea on how many Clusters would be suitable, perform an umap clustering first
residuals_matrix <- assay(combined, "pearson_resid")
residuals_matrix_variable <- residuals_matrix[rownames(residuals_matrix) %in% var_genes, ]
uu.nn <- umap_adt_gex$nn$euclidean$idx
uu.nn[uu.nn == 0] <- 1
## weights as 1-distance
wgt <- 1 - umap_adt_gex$nn$euclidean$dist / min(max(umap_adt_gex$nn$euclidean$dist), 1e5)
wgt[wgt < -1] <- 0
## convert to adjacancy matrix
adj <- matrix(0, ncol(residuals_matrix_variable), ncol(residuals_matrix_variable))
rownames(adj) <- colnames(residuals_matrix_variable)
colnames(adj) <- colnames(residuals_matrix_variable)
for (i in seq_len(ncol(residuals_matrix_variable))) {
  adj[i, colnames(residuals_matrix_variable)[uu.nn[i, ]]] <- wgt[i, ]
}
# convert to weighted graph
g <- graph.adjacency(adj, mode = "undirected", weighted = T, diag = F)
km <- igraph::cluster_louvain(g)
ph.membership <- km$membership
names(ph.membership) <- km$names
umap_cl <- km$membership

# Since UMAP seems to underestimate the number of clusters, add 2 extra
NumberOfClusters <- max(umap_cl) + 2

# Prepare data for bremsc
RNA_count_data <- assay(combined, "counts")
# reduce RNA count data to number of most variable genes specified as command line argument
RNA_count_data <- RNA_count_data[rownames(RNA_count_data) %in% var_genes, ]
# Clustering taking the number of clusters from the umap_clustering + 2)
cluster <- BREMSC(t(adt_counts), RNA_count_data, K = NumberOfClusters, nChains = opt$threads, nMCMC = 500)
cluster.df <- as.data.frame(cluster$clusterID)
cluster.df$barcodes <- colnames(t(adt_counts))
resort_clusterIDs <- as.data.frame(sort(table(cluster$clusterID), decreasing = TRUE))
colnames(resort_clusterIDs) <- c("oldID", "freq")
resort_clusterIDs$newIDs <- rownames(resort_clusterIDs)
cluster.df$newIDs <- resort_clusterIDs$newIDs[match(unlist(cluster$clusterID), resort_clusterIDs$oldID)]

# Plot Clustering
cluster_attributes <- plyr::join(colData_df, cluster.df)
cluster_attributes$cluster <- as.integer(cluster$clusterID)
cluster_attributes$newIDs <- as.integer(cluster_attributes$newIDs)
cluster_umap_final <- ggplot(cluster_attributes, aes(x = ADT_GEX_umap1, y = ADT_GEX_umap2, color = as.factor(newIDs))) +
  geom_point(size = 1) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 15)) +
  scale_color_manual(name = "BREMSC", values = cols33) +
  theme(aspect.ratio = 1) +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 10),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.01, "line"),
    legend.spacing.y = unit(0.1, "line"),
    legend.spacing.x = unit(0.1, "line")
  )
filename <- paste0(opt$outdir, opt$sampleName, "_ADT_GEX__Clustering.png")
ggsave(cluster_umap_final, filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

# generate antibody/expression plots as it is done for the single samples
cat("\n\nCreate antibody / gene expression - plots\n")
gene.list.all <- unique(lookup$Gene)
gene.list.all <- data.frame(Gene = unique(lookup$Gene))

# have matrix with GEX normcounts
# remove all genes that have 0 counts in all cells
normcounts_all.zero.removed <- normcounts(combined)
mask_all_zero <- rowSums(normcounts_all.zero.removed != 0) > 0
normcounts_all.zero.removed <- normcounts_all.zero.removed[mask_all_zero, ]
stopifnot(length(normcounts_all.zero.removed[, 1]) == sum(rowSums(normcounts_all.zero.removed != 0) > 0))

# keep only genes from lookup table that are contained in GEX normcounts matrix
gene.list <- gene.list.all[gene.list.all$Gene %in% rownames(normcounts_all.zero.removed), ]

# sort gene names alphabetically for better overview in plot
gene.list <- sort(gene.list)
gene.list.all <- sort(gene.list.all$Gene)

# get indices of genes in matrix
list_indices_genes <- lapply(seq_along(gene.list), function(x) {
  indices <- match(gene.list[[x]], rownames(normcounts_all.zero.removed))
  indices
})

# Define some general used graphical parameters
fontsize <- theme(
  axis.text = element_text(size = 10),
  axis.title = element_text(size = 10)
)
theme_set(theme_bw(4) + fontsize)

# plot UMAP (based on type_exp), colours = final cell type, as reference for expression plots
p_final_ref <-
  ggplot(colData_df, aes(x = colData_df[[umap1_exp]], y = colData_df[[umap2_exp]], color = celltype_final)) +
  scale_color_manual(
    name = "Cell type",
    values = ct.color[id.final.ct],
    drop = F
  ) +
  geom_point(size = 1) +
  theme(aspect.ratio = 1) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(legend.position = "none")
# plot legend for final cell type UMAP separately
p_legend_ref <-
  ggplot(colData_df, aes(x = colData_df[[umap1_exp]], y = colData_df[[umap2_exp]], color = celltype_final)) +
  scale_color_manual(
    name = "Cell type",
    values = ct.color[id.final.ct],
    drop = F
  ) +
  geom_point(size = 1) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 20)) +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 10),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.01, "line"),
    legend.spacing.y = unit(0.1, "line"),
    legend.spacing.x = unit(0.1, "line")
  )

legend_ref <- cowplot::get_legend(p_legend_ref)

# get all genes that have zero counts in all cells
mask.not.expressed <-
  !(gene.list.all %in% rownames(normcounts_all.zero.removed[unlist(list_indices_genes), , drop = F]))
genes.not.expressed <- sort(gene.list.all[mask.not.expressed])

# get matrix containing only genes of interest
matrix_group <- normcounts_all.zero.removed[unlist(list_indices_genes), , drop = F]

# check first if any genes of interest are expressed in this sample
lookup <- lookup[stringr::str_order(lookup$Antibody, numeric = TRUE), ]
lookup_expressed <- lookup[-which(lookup$Gene %in% genes.not.expressed), ]
lookup_not_expressed <- lookup[which(lookup$Gene %in% genes.not.expressed), ]
rownames(lookup_expressed) <- c(seq_len(length(lookup_expressed$Gene)))
rownames(lookup_not_expressed) <- c(seq_len(length(lookup_not_expressed$Gene)))

# get empty list that will be filled with expression plots per gene
plot_gene <- list()


# iterate over all genes from the lookup table and generate plot with gene
# and antibody expression
numberOfGenes <- length(lookup$Gene)
print("Start iteration:")
for (gene in seq(numberOfGenes)) {
  merged_plots <- list()
  gene_name <- lookup$Gene[gene]
  cat("\nPlotting", as.character(gene_name), "\n")
  # combine colData_df with expression of current gene
  colData_df_gene <- colData_df

  ###   differentiate between genes with RNA-expression and without
  colour_when_zero <- NULL
  # if gene is expressed add normcounts to colData_df_gene
  if (gene_name %in% rownames(matrix_group)) {
    colData_df_gene$normcounts <- matrix_group[as.character(gene_name), ]
    colour_when_zero <- "slateblue4"
    # if gene is not expressed (not in expression matrix) add zeroes to colData_df_gene
  } else {
    colData_df_gene$normcounts <- rep.int(0, length(colData_df_gene$barcodes))
    colour_when_zero <- "gray"
  }
  # maximum count found for current gene
  max_count <- max(colData_df_gene$normcounts)
  # upper limit of gene expression colour scale is either maximum count or 3
  upper_limit <- max(3, max_count)
  # Plotting both, gene and antibody name, only if they are not identical
  if (toupper(gene_name) == toupper(lookup[gene, ]$Antibody)) {
    legend <- gene_name
  } else {
    legend <- paste(gene_name, " / ", lookup[gene, ]$Antibody, sep = "")
  }

  ###   plot gene expression   ###
  genes_zero_expression <- colData_df_gene[which(colData_df_gene$normcounts == 0), ]
  genes_nonZero_expression <- colData_df_gene[which(colData_df_gene$normcounts > 0), ]
  print("Plot gene expression:")
  expr <-
    ggplot(colData_df_gene, aes(x = colData_df_gene[[umap1_exp]], y = colData_df_gene[[umap2_exp]])) +
    # all cells (obsolete?)
    geom_point(colour = "gray", size = rel(0.001)) +
    # cells with zero expression
    geom_point(
      data = genes_zero_expression,
      x = genes_zero_expression[[umap1_exp]],
      y = genes_zero_expression[[umap2_exp]],
      size = rel(0.001),
      colour = colour_when_zero
    ) +
    # cells with non-zero expression
    geom_point(
      data = genes_nonZero_expression,
      x = genes_nonZero_expression[[umap1_exp]],
      y = genes_nonZero_expression[[umap2_exp]],
      aes(color = normcounts),
      size = rel(0.001)
    ) +
    scale_color_gradientn(
      name = "counts",
      na.value = "gray",
      colours = c("slateblue3", "royalblue1", "aquamarine3", "khaki", 383, "sienna1", "orangered4"),
      limits = c(1, max(3, upper_limit)),
      breaks = c(floor(upper_limit / 3), round(2 * (upper_limit / 3)), upper_limit)
    ) +
    coord_fixed(ratio = 1) +
    ggtitle(legend) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      plot.title = element_text(
        face = "bold", hjust = 0.5, vjust = 0.1, size = 14
      ),
      plot.margin = margin(1, 1, 1, 1, "pt"),
      legend.key.size = unit(0.8, "line"),
      legend.text = element_text(size = rel(3)),
      legend.title = element_text(size = rel(3)),
      legend.margin = margin(0.2, 0.2, 0.2, 0.2),
      axis.text = element_text(size = rel(3)),
      aspect.ratio = 1
    )
  print("Plot antibody expression:")
  # combine colData_df with antibody
  my_ab <- paste0("Norm_ADT_", lookup[gene, ]$Antibody)
  colData_df_antibody <- colData_df
  colData_df_antibody$expressed <- as.numeric(colData_df[, my_ab])
  # take natural (ln) log of ADT counts
  colData_df_antibody$expressed <- round(log(colData_df_antibody$expressed))
  ab <- lookup[gene, ]$Antibody

  # plot antibody plot
  # Define the lowest threshold accross all samples for the antibody to be the reference for gray cells.
  min_threshold <- min(thresholds[lookup[gene, ]$Antibody, ])
  upper_limit <- max(min_threshold * 3, max(colData_df_antibody$expressed), 3, na.rm = TRUE)
  medium_break <- (min_threshold + upper_limit) / 2
  upper_break <- (medium_break + upper_limit) / 2
  lower_break <- (medium_break + min_threshold) / 2
  ###   differentiate between antibodies that are expressed and those that are not
  # if antibody is not expressed
  colData_df_antibody$expressed[is.na(colData_df_antibody$expressed)] <- 0
  if (all(colData_df_antibody$expressed == 0) | all(is.na(colData_df_antibody$expressed))) {
    if (all(is.na(colData_df_antibody$expressed))) {
      colData_df_antibody$expressed <- 0
    }
    antibody_plot <- ggplot(
      colData_df_antibody,
      aes(x = colData_df[[umap1_exp]], y = colData_df[[umap2_exp]], color = expressed)
    ) +
      geom_point(
        data = colData_df_antibody[which(colData_df_antibody$expressed == 0), ],
        size = rel(0.001),
        colour = "gray"
      ) +
      scale_color_gradientn(
        name = "ln counts",
        na.value = "gray",
        colours = c("slateblue4", "royalblue1", "aquamarine3", "khaki", 383, "sienna1", "orangered4"),
        #        limits = c(min(thresholds[lookup[gene, ]$Antibody, ]), upper_limit),
        limits = c(0, upper_limit),
        breaks = c(0, upper_limit / 2, upper_limit),
        # breaks = c(min(thresholds[lookup[gene, ]$Antibody, ]), medium_break, upper_limit),
        labels = c(
          paste0("<", min(thresholds[lookup[gene, ]$Antibody, ])),
          medium_break,
          upper_limit
        )
      ) +
      coord_fixed(ratio = 1) +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(
          face = "bold",
          hjust = 0.5,
          vjust = 0.1,
          size = 8
        ),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "pt"),
        legend.key.size = unit(0.8, "line"),
        legend.text = element_text(size = rel(3)),
        legend.title = element_text(size = rel(3)),
        legend.margin = margin(0.2, 0.2, 0.2, 0.2),
        axis.text = element_text(size = rel(3)),
        aspect.ratio = 1
      ) +
      guides(fill = guide_legend(title = "counts"))
    # if antibody is expressed
  } else {
    antibody_plot <- ggplot(
      colData_df_antibody,
      aes(x = colData_df[[umap1_exp]], y = colData_df[[umap2_exp]], color = expressed)
    ) +
      # geom_point(data = colData_df_antibody[colData_df_antibody$expressed == 0, ],
      #          size = rel(0.001),
      #         colour = "gray") +
      # geom_point(data = colData_df_antibody[which(colData_df_antibody$expressed > thresholds[ab, ]), ],
      geom_point(
        data = colData_df_antibody,
        #     na.rm = TRUE,
        size = rel(0.001)
      ) +
      # geom_point(data = colData_df_antibody[which(colData_df_antibody$expressed > max(thresholds[ab, ], upper_limit)), ],
      #           size = rel(0.001),
      #           colour = "orangered4") +
      scale_color_gradientn(
        name = "ln counts",
        na.value = "gray",
        colours = c("slateblue4", "royalblue1", "aquamarine3", "khaki", 383, "sienna1", "orangered4"),
        # limits = c(min(thresholds[lookup[gene, ]$Antibody, ]), upper_limit),
        # breaks = c(min(thresholds[lookup[gene, ]$Antibody, ]), medium_break, upper_limit),
        limits = c(0, upper_limit),
        breaks = c(0, upper_limit / 2, upper_limit),
        labels = c(
          paste0("<", min(thresholds[lookup[gene, ]$Antibody, ])),
          medium_break,
          upper_limit
        )
      ) +
      coord_fixed(ratio = 1) +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(
          face = "bold", hjust = 0.5, vjust = 0.1, size = 8
        ),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "pt"),
        legend.key.size = unit(0.8, "line"),
        legend.text = element_text(size = rel(3)),
        legend.title = element_text(size = rel(3)),
        legend.margin = margin(0.2, 0.2, 0.2, 0.2),
        axis.text = element_text(size = rel(3)),
        aspect.ratio = 1
      ) +
      guides(fill = guide_legend(title = "counts"))
  }
  # combine expression plot and antibody plot of current gene to one plot
  plot_gene[[gene]] <- expr / antibody_plot
}

# TODO What does this syntax mean?
layout <- c(
  area(1, 2, 5, 6),
  area(1, 7, 5, 11),
  area(1, 12, 5, 16),
  area(1, 16, 5, 19),
  area(6, 1, 8, 4),
  area(6, 4, 8, 7),
  area(6, 7, 8, 10),
  area(6, 10, 8, 13),
  area(6, 13, 8, 16),
  area(6, 16, 8, 19),
  area(9, 1, 11, 4),
  area(9, 4, 11, 7),
  area(9, 7, 11, 10),
  area(9, 10, 11, 13),
  area(9, 13, 11, 16),
  area(9, 16, 11, 19)
)

###   Save expression plots   ###
cat("\nSaving Expression Plots:\n")
numberOfPlots <- length(plot_gene)
# only plot 12 plots into one file
numberOfPages <- numberOfPlots %/% 12
if (numberOfPlots %% 12 != 0) {
  numberOfPages <- numberOfPages + 1
}

dir.create(paste0(opt$outdir, "ExpressionPlots/"))
print("Number of subplots to plot")
print(numberOfPlots)

# save plots
for (i in seq(1, numberOfPlots, 12)) {
  plot_collection <- p_final_ref + legend_ref + umap_sample + legend_sample
  for (p in seq(1, 12)) {
    plot_collection <- plot_collection + plot_gene[p + i]
  }
  plot_collection <- plot_collection + plot_layout(design = layout)
  final <- plot_collection + plot_annotation(caption = paste("Page ", (i + 12 - 1) / 12), " of ", numberOfPages) & theme(plot.caption = element_text(size = 10))
  out_expression_plots <- paste0(opt$outdir, "ExpressionPlots/", opt$sampleName, "_", type_exp, "_")
  if (i + 12 > numberOfPlots) {
    ggplot2::ggsave(
      filename = paste(out_expression_plots, lookup$Antibody[i + 1], "to", lookup$Antibody[numberOfPlots], "plot.png", sep = "_"),
      width = 21,
      height = 14,
      plot = final
    )
  } else if (i + 12 == numberOfPlots) {
    ggplot2::ggsave(
      filename = paste(out_expression_plots, lookup$Antibody[i + 1], "to", lookup$Antibody[numberOfPlots], "plot.png", sep = "_"),
      width = 21,
      height = 14,
      plot = final
    )
    # all plots have been plotted already
    break
  } else {
    ggplot2::ggsave(
      filename = paste(out_expression_plots, lookup$Antibody[i + 1], "to", lookup$Antibody[i + 12], "plot.png", sep = "_"),
      width = 21,
      height = 14,
      plot = final
    )
  }
}

cat("\n\n\nPrint sessionInfo():\n\n")
print(sessionInfo())
