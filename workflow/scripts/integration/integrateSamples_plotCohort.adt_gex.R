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
   #make_option("--cohort_list", type = "character", help = "List with paths to the input RDS files that contain SCE objects with GEX results. INPORTANT: GEX and GSVA input files have to be in the same order."),
   make_option("--cohort_object", type = "character", help = "Integrated Seurat Object"),
   #make_option("--gsva_list", type = "character", help = "List with paths to the input RDS files that contain SCE objects with GSVA results. IMPORTANT: GEX and GSVA input files have to be in the same order."),
   make_option("--colour_config", type = "character", help = "Path to config file table that sorts colours to cell types"),
   #make_option("--selectedGenes", type = "character", help = "Tab delimited text file with the genes of interest in the first column. No header. The genes will be plotted and given out in chunks, alphabetical order."),
   make_option("--add_meta", type = "character", help = "Tab delimited text file with additional meta data. The first column must have header 'sampleID'. For a default run the full file name of the respective SCE analysis object is required, refer to script documentation on different sampleID formatting. Column headers describe categories. One line per input sample is required. Meta data file itself is optional: without table simply do not use --add_meta"),
   make_option("--lookup", type = "character", help = "Full path to a lookup table matching protein to antibody names."),
   make_option("--thresholds", type = "character", help = "Full path to a table listing adt thresholds for every sample included in the annalysis."),
   make_option("--outdir", type = "character", help = "Full path to output directory."),
   make_option("--outName", type = "character", help = "Prefix name of output files.")
 )
 opt_parser <- OptionParser(option_list = option_list)
 opt <- parse_args(opt_parser)


# convenience function for string concatenation
"%&%" <- function(a, b) paste(a, b, sep = "")
empty_string <- ""


###################################
###   Set plotting Parameters   ###
###################################

# If one of the values plot_aspect_ratio, chunksize, or number of columns is changed in the following section
# the other values should be checked as well.
# Also, check physical size of the generated output file in the respective plotting section of the script.

# IMPORTANT: when generating UMAP plots with the function RunUMAP() from Seurat the parameters "spread", "min.dist" and "n.neighbors"
# have a great impact on the resulting plot. They should be checked and optimised for each project.

# for all UMAP plots the same aspect ratio is used to have a comparable visualisation
plot_aspect_ratio <- 0.7

### The selected genes and gene sets are plotted in chunks so that the file size is not too large:

# chunksize_genes determines the number of genes shown in one expression plot
chunksize_genes <- 9
# columns_genes determines the number of columns the gene expression plots are arranged in.
## E.g. "columns_genes <- 1" arranges all plots in the output file below one another in a single column.
columns_genes <- 3

# chunk size of gene sets in GSVA plots
chunksize_genesets <- 12
# number of columns in GSVA plots
columns_genesets <- 4

# override the theme change of cowplot
theme_set(theme_grey())
print(Sys.time())


########################
###   Read in data   ###
########################

# read in color config
#config <- read.csv("colour_config.merged.txt", sep = "\t", stringsAsFactors = FALSE)
config <- read.csv(opt$colour_config, sep = "\t", stringsAsFactors = FALSE)
print(config)

# read in seurat object
#seurat_integrated <- readRDS("Cohort.integrated_seurat_object.RDS")
seurat_integrated <- readRDS(opt$cohort_object)
print(seurat_integrated)
cell_attributes <- as.data.frame(Embeddings(seurat_integrated, reduction = "wnn.umap"))
print(cell_attributes)
print(rownames(cell_attributes))
print(substring(rownames(cell_attributes),0,16))
print(unique(substring(rownames(cell_attributes),0,16)))
print(length(rownames(cell_attributes)))
print(length(unique(substring(rownames(cell_attributes),0,16))))
print(duplicated(substring(rownames(cell_attributes),0,16)))
cell_attributes <- cell_attributes[!duplicated(substring(rownames(cell_attributes),0,16)), ]
print(cell_attributes)
rownames(cell_attributes) <- substring(rownames(cell_attributes),0,16)
print(head(cell_attributes))


metadata <- seurat_integrated@meta.data
metadata <- metadata[!duplicated(substring(rownames(metadata),0,16)),]
print(head(rownames(metadata)))
rownames(metadata) <- substring(rownames(metadata),0,16)
print(head(rownames(metadata)))
print(colnames(metadata))
metadata$barcodes <- NULL
metadata$barcodes <- rownames(metadata)
print(colnames(metadata))
print(metadata$barcodes)
print(head(metadata))
cell_attributes$barcodes <- NULL
cell_attributes$barcodes <- rownames(cell_attributes)
cell_attributes <- plyr::join(cell_attributes, metadata, by="barcodes")
rownames(cell_attributes) <- cell_attributes$barcodes
cell_attributes$Row.names <- NULL
print(head(cell_attributes))



print(head(metadata))
number_samples <- length(unique(metadata$sampleID))
# optionally, read in an additional meta data table
print(empty_string)
 if (is.character(opt$add_meta)) {
#   addMeta_table <- read.csv("metadata_table.PBMC_Figure2", sep = "\t", stringsAsFactors = FALSE)
   addMeta_table <- read.csv(opt$add_meta, sep = "\t", stringsAsFactors = FALSE)
   ## Use sampleID up to (but excluding) the second "." as shorter sampleID
   addMeta_table$sampleID <- regmatches(addMeta_table$sampleID, regexpr("[^\\.]+[^\\.]+", addMeta_table$sampleID))
   print("Reading in meta data table: check if input samples and meta data sampleIDs are identical")
   print(addMeta_table$sampleID)
   stopifnot(sort(addMeta_table$sampleID) == sort(unique(metadata$sampleID)))
 } else if (is.null(opt$add_meta)) {
   print("No input file with additional meta data given.")
   addMeta_table <- NULL
 } else {
   print("Something went wrong reading in the meta data table.")
 }

# Read in Marker Lists
# if (is.character(opt$markerGenes)) {
# MarkerGenes <-  read.csv(opt$markerGenes, sep = "\t", stringsAsFactors = FALSE)
# } else if (is.null(opt$markerGenes)) {
#   print("No input file with additional meta data given.")
#   markerGenes <- NULL
# }
#print(markerGenes)
#if (is.character(opt$markerADTs)) {
# MarkerADTs <-  read.csv(opt$markerADTs, sep = "\t", stringsAsFactors = FALSE)
# } else if (is.null(opt$markerADTs)) {
#   print("No input file with additional meta data given.")
#   markerADTs <- NULL
# }
#print(markerADTs)

# Read in lookup table

lookup <-read.csv(opt$lookup, sep = "\t", stringsAsFactors = FALSE)
lookup$Antibody <- sub("_", "-", lookup$Antibody)
print(lookup)
# Read in thresholds
thresholds <- read.csv(opt$thresholds, sep= "\t", stringsAsFactors = FALSE)
thresholds$Antibody<- sub("_", "-", thresholds$Antibody)
rownames(thresholds) <- toupper(thresholds$Antibody)
print(tail(thresholds,20))

set.seed(2)

config <- read.csv(opt$colour_config, sep = "\t", stringsAsFactors = FALSE)

config$class <- NULL
config$indication <- NULL
# have character vector of colour names from the colour config

config <- rbind(config, c("uncertain", "grey50", "uncertain"))
config <- rbind(config, c("unknown", "black", "unknown"))
# give names to colours from other column of colour config
ct.color <- config$colour
names(ct.color) <- config$mapping
print(ct.color)
# This is done in the integration script currently. Think about were it makes more sense
# make sure 'celltype_major' and 'celltype_final' metadata columns contain factors and are mapped to the correct celltype
#seurat_integrated[["celltype_major"]] <- as.factor(config$mapping[match(as.factor(seurat_integrated[[]]$celltype_major), config$cell_type)])
#seurat_integrated[["celltype_final"]] <- as.factor(config$mapping[match(as.factor(seurat_integrated[[]]$celltype_final), config$cell_type)])

# make sure "original objects" contain the correct mapping for the celltypes

# get integer vector containing the index of the matching colour in ct.color for each
# level of seurat_integrated[[]]$celltype_major
id.first.ct <- match(levels(seurat_integrated[[]]$celltype_major), names(ct.color))
id.first.ct
id.final.ct <- match(levels(seurat_integrated[[]]$celltype_final), names(ct.color))
id.final.ct

colstypes <- c("red2", "green4", "green2", "mediumorchid2", "blue2", "cyan2", "yellow1", "purple", "brown",
               "chocolate1", "chartreuse2", "darkgoldenrod3", "steelblue1", "slateblue3", "olivedrab4", "gold2",
               "violetred3", "darkcyan", "orchid3", "darksalmon", "darkslategrey", "khaki", "indianred2", "magenta", "slategray2",
               "olivedrab1", "mediumaquamarine", "hotpink", "yellow3",
               "bisque4", "darkseagreen1", "dodgerblue3",
               "deeppink4", "sienna4", "mediumorchid4")
###   Plot general cohort composition   ###
print(empty_string)
print("Plot general cohort composition.")
print(Sys.time())


p1 <- DimPlot(seurat_integrated, reduction = 'wnn.umap', , group.by = "sample_name", repel = TRUE) +
  theme(aspect.ratio = plot_aspect_ratio)
p1 

ggsave(filename = paste0(opt$outdir, opt$outName, ".sample_integration_all_samples.png"),
       width = 28, height = (22 + number_samples*0.3), dpi = 300, units = "cm")

p1.gex <- DimPlot(seurat_integrated, reduction = 'rna.umap', group.by = 'sample_name',  
              repel = TRUE) + NoLegend() + theme(aspect.ratio = plot_aspect_ratio)
p1.adt <- DimPlot(seurat_integrated, reduction = 'adt.umap', group.by = 'sample_name', 
              repel = TRUE) + theme(aspect.ratio = plot_aspect_ratio)
p1.gex + p1.adt
ggsave(filename = paste0(opt$outdir, opt$outName, ".sample_integration_all_samples.individual.png"),
       width = 28, height = (22 + number_samples*0.3), dpi = 300, units = "cm")

###
p2 <- DimPlot(seurat_integrated, reduction = 'wnn.umap', label = FALSE,group.by = "celltype_final",cols = ct.color[id.final.ct],  repel = TRUE, label.size = 2.5) +
  theme(aspect.ratio = plot_aspect_ratio) +
  guides(color = guide_legend(title = "celltype_final", override.aes = list(size = 3)))
p2 

ggsave(filename = paste0(opt$outdir, opt$outName, ".sample_integration_celltypes_final.png"),
       width = 34, height = 20, dpi = 300, units = "cm")
p2.gex <- DimPlot(seurat_integrated, reduction = 'rna.umap', group.by = 'celltype_final',  
                  repel = TRUE) + NoLegend() + theme(aspect.ratio = plot_aspect_ratio)
p2.adt <- DimPlot(seurat_integrated, reduction = 'adt.umap', group.by = 'celltype_final', 
                  repel = TRUE) + theme(aspect.ratio = plot_aspect_ratio)
p2.gex + p2.adt
ggsave(filename = paste0(opt$outdir, opt$outName, ".sample_integration_all_samples_celltypes_final.individual.png"),
       width = 28, height = (22 + number_samples*0.3), dpi = 300, units = "cm")

print("Written plot 2")
###
p3 <- DimPlot(seurat_integrated, reduction = 'wnn.umap', label = FALSE,group.by = "celltype_major",cols = ct.color[id.final.ct],  repel = TRUE, label.size = 2.5) +
  theme(aspect.ratio = plot_aspect_ratio) +
  guides(color = guide_legend(title = "celltype_major", override.aes = list(size = 3)))
p3 

ggsave(filename = paste0(opt$outdir, opt$outName, ".sample_integration_celltypes_major.png"),
       width = 34, height = 20, dpi = 300, units = "cm")
p3.gex <- DimPlot(seurat_integrated, reduction = 'rna.umap', group.by = 'celltype_major',cols = ct.color[id.final.ct],  
                  repel = TRUE) + NoLegend() + theme(aspect.ratio = plot_aspect_ratio)
p3.adt <- DimPlot(seurat_integrated, reduction = 'adt.umap', group.by = 'celltype_major', cols = ct.color[id.final.ct],
                  repel = TRUE) + theme(aspect.ratio = plot_aspect_ratio)
p3.gex + p3.adt
ggsave(filename = paste0(opt$outdir, opt$outName, ".sample_integration_all_samples_major.individual.png"),
       width = 28, height = (22 + number_samples*0.3), dpi = 300, units = "cm")

print("Written plot 3")
p4 <- DimPlot(seurat_integrated, reduction = "wnn.umap", group.by = "seurat_clusters", cols = colstypes, label = TRUE) +
  theme(aspect.ratio = plot_aspect_ratio) +
  guides(color = guide_legend(title = "seurat_clusters", override.aes = list(size = 3)))
p4
ggsave(filename = paste0(opt$outdir, opt$outName, ".sample_integration_seurat_umap_clusters.png"),
       width = 30, height = 20, dpi = 300, units = "cm")
print("Written plot 4")


###   Characterise cohort   ###
print(empty_string)
print("Characterise cohort plots.")
# Plot cycle phase
p22 <- DimPlot(seurat_integrated, reduction = "wnn.umap", group.by = "cycle_phase", cols = colstypes) +
  theme(aspect.ratio = plot_aspect_ratio) +
  guides(color = guide_legend(title = "cycle_phase", override.aes = list(size = 3)))
p22
ggsave(filename = paste0(opt$outdir, opt$outName, ".cellcycle_phase.png"), width = 30, height = 20, dpi = 300, units = "cm")

# Plot fraction of MT reads
p23 <- FeaturePlot(seurat_integrated, reduction = "wnn.umap", features = c("fractionMT")) +
  theme(aspect.ratio = plot_aspect_ratio)
p23
ggsave(filename = paste0(opt$outdir, opt$outName, ".fraction_MT.png"), width = 30, height = 20, dpi = 300, units = "cm")

# Plot n_umi
p24 <- FeaturePlot(seurat_integrated,reduction = "wnn.umap", features = c("n_umi")) +
  theme(aspect.ratio = plot_aspect_ratio)
p24
ggsave(filename = paste0(opt$outdir, opt$outName, ".n_umi.png"), width = 30, height = 20, dpi = 300, units = "cm")

# Plot n_gene
p25 <- FeaturePlot(seurat_integrated,reduction = "wnn.umap", features = c("n_gene")) +
  theme(aspect.ratio = plot_aspect_ratio)
p25
ggsave(filename = paste0(opt$outdir, opt$outName, ".n_gene.png"), width = 30, height = 20, dpi = 300, units = "cm")

# Plot log_umi
p26 <- FeaturePlot(seurat_integrated, reduction = "wnn.umap",features = c("log_umi")) +
  theme(aspect.ratio = plot_aspect_ratio)
p26
ggsave(filename = paste0(opt$outdir, opt$outName, ".log_umi.png"), width = 30, height = 20, dpi = 300, units = "cm")

# Plot g2m_score
p27 <- FeaturePlot(seurat_integrated,reduction = "wnn.umap", features = c("g2m_score")) +
  theme(aspect.ratio = plot_aspect_ratio)
ggsave(filename = paste0(opt$outdir, opt$outName, ".g2m_score.png"), width = 30, height = 20, dpi = 300, units = "cm")

# Plot s_score
p28 <- FeaturePlot(seurat_integrated,reduction = "wnn.umap", features = c("s_score")) +
  theme(aspect.ratio = plot_aspect_ratio)
p28
ggsave(filename = paste0(opt$outdir, opt$outName, ".s_score.png"), width = 30, height = 20, dpi = 300, units = "cm")

###############
# Plot Metadata
###############
print(" Plot Metadata")
pm1 <- DimPlot(seurat_integrated,reduction = "wnn.umap", group.by  = c("hashtag")) +
  theme(aspect.ratio = plot_aspect_ratio) +
  guides(color = guide_legend(title = "hashtag", override.aes = list(size = 3)))
pm1
ggsave(filename = paste0(opt$outdir, opt$outName, ".hashtag.png"), width = 30, height = 20, dpi = 300, units = "cm")



pm2 <- DimPlot(seurat_integrated,reduction = "wnn.umap", group.by  = c("dissociation")) +
  theme(aspect.ratio = plot_aspect_ratio) +
  guides(color = guide_legend(title = "Digesting method", override.aes = list(size = 3)))
pm2
ggsave(filename = paste0(opt$outdir, opt$outName, ".digestion.png"), width = 30, height = 20, dpi = 300, units = "cm")


#######################################
# Plot AB Expression
#######################################


# generate antibody/expression plots as it is done for the single samples
cat("\n\nCreate antibody / gene expression - plots\n")

  
my_merge <- function(df1, df2){ 
  merge(df1, df2, by = "rownames", all=TRUE)
}

#joined_gene_expr <- matrix_list %>%  purrr::reduce(my_merge)
#rownames(joined_gene_expr) <- joined_gene_expr$rownames
#joined_gene_expr$rownames <- NULL
#joined_gene_expr <- as.data.frame(t(joined_gene_expr))


DefaultAssay(seurat_integrated) <- "adt"

print(rownames(cell_attributes))

cell_attributes <- as.data.frame(Embeddings(seurat_integrated, reduction = "wnn.umap"))
cell_attributes$barcodes<- rownames(cell_attributes)
cell_attributes <- merge(cell_attributes, as.data.frame(seurat_integrated@meta.data), by=0)  
rownames(cell_attributes) <- cell_attributes$barcodes
cell_attributes$Row.names <- NULL

p_final_ref <-
  ggplot(cell_attributes, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = celltype_final)) +
  scale_color_manual(name = "Cell type",
                     values = ct.color[id.final.ct],
                     drop = F) +
  geom_point(size = 1) +
  theme(aspect.ratio = 1,panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("wnnUMAP 1") +
  ylab("wnnUMAP 2") +
  theme(legend.position = "none")
p_final_ref
p_legend_ref <-
  ggplot(cell_attributes, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = celltype_final)) +
  scale_color_manual(name = "Cell type",
                     values = ct.color[id.final.ct],
                     drop = F) +
  geom_point(size = 1) +
  xlab("wnnUMAP 1") +
  ylab("wnnUMAP 2") +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 15)) +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 10),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.01, "line"),
    legend.spacing.y = unit(0.1, "line"),
    legend.spacing.x = unit(0.1, "line"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
  )

p_legend_ref 

legend_ref <- cowplot::get_legend(p_legend_ref)

plot_gene <- list()
# iterate over all genes from the lookup table and generate plot with gene
# and antibody expression
numberOfGenes <- length(lookup$Gene)


for (gene in seq(numberOfGenes)) {
  print(gene)
  merged_plots <- list()
  gene_name <- lookup$Gene[gene]
  cat("\nPlotting", as.character(gene_name), "\n")
  # Plotting both, gene and antibody name, only if they are not identical
  if (toupper(gene_name) == toupper(lookup[gene, ]$Antibody)) {
    legend <- gene_name
  } else {
    legend <- paste(gene_name, " / ", lookup[gene, ]$Antibody, sep = "")
  }
  if (!(gene_name %in% rownames(seurat_integrated@assays$SCT@counts))){
    print(paste0(gene_name," not detected."))
    expr <- ggplot(cell_attributes, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
      # all cells (obsolete?)
      geom_point(color="gray", size = 0.1)+
      ggtitle(legend) +theme_classic() +
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
          size = 16
        ),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        legend.key.size = unit(0.8, "line"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.margin = margin(0.2, 0.2, 0.2, 0.2),
        axis.text = element_text(size = 5),
        aspect.ratio = 1
      ) 
  } else {
    print(paste0(gene_name," detected."))
    ## combine colData_df with expression of current gene
   # colData_df_gene <- cell_attributes[order(cell_attributes$barcodes),]
  #  ordered_gene_expr <- as.data.frame(joined_gene_expr[order(rownames(joined_gene_expr)),])
  #  colnames(ordered_gene_expr) <- toupper(colnames(ordered_gene_expr))
  #  colData_df_gene$expressed <- ordered_gene_expr[,gene_name]
    ###   differentiate between genes with RNA-expression and without
    colour_when_zero <- NULL
    # if gene is expressed add normcounts to colData_df_gene
    if (max(seurat_integrated@assays$SCT@counts[gene_name,], na.rm = TRUE)>0) {
      colour_when_zero <- "slateblue4"
      # if gene is not expressed (not in expression matrix) add zeroes to colData_df_gene
    } else {
      colour_when_zero <- "gray"
    }
    print(colour_when_zero)
    # maximum count found for current gene
    max_count <- max(seurat_integrated@assays$SCT@counts[gene_name,])
    # upper limit of gene expression colour scale is either maximum count or 3
    upper_limit <- max(3, max_count)
  ###   plot gene expression   ###
  #genes_zero_expression <- colData_df_gene[which(colData_df_gene$expressed == 0), ]
   genes_zero_expression <- as.data.frame(seurat_integrated@assays$SCT@counts[gene_name, which(seurat_integrated@assays$SCT@counts[gene_name,]==0), ])
  my_rownames <- rownames(genes_zero_expression)
  my_new_rownames <- substring(my_rownames[!duplicated(substring(rownames(genes_zero_expression),0,16))], 0,16)
  colnames(genes_zero_expression) <- "expressed"
  genes_zero_expression <- as.data.frame(genes_zero_expression[!duplicated(substring(rownames(genes_zero_expression),0,16)), ])
  colnames(genes_zero_expression) <- "expressed"
  rownames(genes_zero_expression) <- my_new_rownames
  colnames(genes_zero_expression) <- "expressed"
  genes_zero_expression <- merge(cell_attributes, genes_zero_expression, by=0)
  genes_nonZero_expression <- as.data.frame(seurat_integrated@assays$SCT@counts[gene_name, which(seurat_integrated@assays$SCT@counts[gene_name,]> 0), ])
  my_rownames <- rownames(genes_nonZero_expression)
  my_new_rownames <- substring(my_rownames[!duplicated(substring(rownames(genes_nonZero_expression),0,16))], 0,16)
  genes_nonZero_expression <- as.data.frame(genes_nonZero_expression[!duplicated(substring(rownames(genes_nonZero_expression),0,16)), ])

  rownames(genes_nonZero_expression) <- my_new_rownames
  colnames(genes_nonZero_expression) <- "expressed"
  genes_nonZero_expression <- merge(cell_attributes, genes_nonZero_expression, by=0)
   DefaultAssay(seurat_integrated) <- "SCT"
  expr <- FeaturePlot(seurat_integrated, features = gene_name, slot="counts", 
                     reduction = 'wnn.umap', max.cutoff = "q99",  ncol = 1)
  expr
  expr <- expr + 
    geom_point(data = genes_zero_expression,
               aes(x = wnnUMAP_1,
                   y = wnnUMAP_2),
               colour = "slateblue4", size = 0.1) +
    # cells with non-zero expression
    geom_point(data = genes_nonZero_expression,
               aes(x = wnnUMAP_1,
                   y = wnnUMAP_2,
                   color = expressed,
               ),size = 0.1)+
    scale_color_gradientn(
      name = "counts",
      na.value = "gray",
      colours = c("slateblue3", "royalblue1", "aquamarine3", "khaki", 383, "sienna1", "orangered4"),
      limits = c(1, max(3, upper_limit)),
      breaks = c(floor(upper_limit / 3), round(2 * (upper_limit / 3)), upper_limit)
    ) +
    ggtitle(legend) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      #   face = "bold", hjust = 0.5, vjust = 0.1, size = 14),
      plot.margin = margin(1, 1, 1, 1, "pt"),
      legend.key.size = unit(0.8, "line"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.margin = margin(0.2, 0.2, 0.2, 0.2),
      axis.text = element_text(size = 5),
      aspect.ratio = 1,
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        vjust = 0.1,
        size = 16
      ),
    )
  expr
  }
  # combine colData_df with antibody
  my_ab <- lookup[gene, ]$Antibody
  my_ab_df <- t(as.data.frame(seurat_integrated@assays$adt@counts))
  colnames(my_ab_df) <- toupper(colnames(my_ab_df))
  rownames(cell_attributes) <- cell_attributes$barcodes.x
  colData_df_antibody <- merge(cell_attributes, my_ab_df, by=0)  
  rownames(colData_df_antibody) <- colData_df_antibody$Row.names
  if (paste0("ADT-",toupper(my_ab)) %in% colnames(my_ab_df)){
    colData_df_antibody$expressed <- log(colData_df_antibody[,paste0("ADT-",toupper(my_ab))]+1)
  }  else {
    colData_df_antibody$expressed <- 0}
  print(head(my_ab_df[,paste0("ADT-",toupper(my_ab))]))
  print(head(my_ab_df[,paste0("ADT-",toupper(my_ab))]+1))
  print(head(colData_df_antibody$expressed))
  #
  colData_df_antibody$expressed_thresholds <- colData_df_antibody$expressed
  min_threshold <- 0
  for (sample in unique(colData_df_antibody$sampleID)){
    my_threshold <- thresholds[toupper(my_ab),sample]
    print(my_threshold)
    if (my_threshold >0 & my_threshold < min_threshold){
      min_threshold <- my_threshold
      print("min_threshold")
      print(min_threshold)
    }
    set_to_na <- which(colData_df_antibody[colData_df_antibody$sampleID==sample,]$expressed<my_threshold)
    if (length(set_to_na) >0){
    colData_df_antibody[colData_df_antibody$sampleID==sample,]$expressed_thresholds[set_to_na] <- NA
    }
  }
  colData_df_antibody$expressed <- colData_df_antibody$expressed_thresholds
  print(head(colData_df_antibody$expressed_thresholds))
  #if (colData_df_antibody$expressed < thresholds[])
# take natural (ln) log of ADT counts
  # Define the lowest threshold accross all samples for the antibody to be the reference for gray cells.
  #temp_threshold <- thresholds
  #temp_threshold$Antibody <- NULL
  #rownames(temp_threshold) <- toupper(rownames(temp_threshold))
  #min_threshold <- min(temp_threshold[toupper(my_ab),unique(colData_df_antibody$sampleID)])
  upper_limit <- round(max(min_threshold * 3, max(colData_df_antibody$expressed, na.rm=TRUE), 3, na.rm = TRUE),1)
  medium_break <- round((min_threshold + upper_limit) / 2,1)
 # upper_break <- round((medium_break + upper_limit) / 2,1)
  lower_break <- round((medium_break + min_threshold) / 2,1)
#  colData_df_antibody$expressed <- round(log(colData_df_antibody$expressed))
  
  ###   differentiate between antibodies that are expressed and those that are not
  # if antibody is not expressed
  if (all(colData_df_antibody$expressed == 0)) {
    print(paste0(my_ab," ab not detected."))
    print(head(colData_df_antibody))
    antibody_plot <- ggplot(as.data.frame(colData_df_antibody),
                            aes(x = wnnUMAP_1, y = wnnUMAP_2, color = expressed)) +
      geom_point(data = colData_df_antibody$expressed,
                 size = 0.1,
                 colour = "gray") +
      #scale_color_gradientn(
      #  name = "ln counts",
    #    na.value = "slateblue4",
     #   colours = c("slateblue3", "royalblue1", "aquamarine3", "khaki", 383, "sienna1", "orangered4"),
      #  limits = c(0, upper_limit),
       # breaks = c(0, upper_limit / 2, upper_limit),
        # breaks = c(min(thresholds[lookup[gene, ]$Antibody, ]), medium_break, upper_limit),
      # labels = c(paste0("<", min_threshold),
       #   medium_break,
        #  upper_limit)
    #  ) +
      coord_fixed(ratio = 1) +
      theme_classic() +
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
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.margin = margin(0.2, 0.2, 0.2, 0.2),
        axis.text = element_text(size = 5),
        aspect.ratio = 1
      ) +
      guides(fill = guide_legend(title = "counts"))
    # if antibody is expressed
  } else {
    DefaultAssay(seurat_integrated) <- "adt"
    print(rownames(seurat_integrated@assays$adt@counts))
    antibody <- FeaturePlot(seurat_integrated, features = paste0("ADT-",my_ab), slot="counts", 
                            reduction = 'wnn.umap', max.cutoff = "q99",  ncol = 1)
    print(antibody)
    antibody <- antibody + 
      geom_point(data = colData_df_antibody[is.na(colData_df_antibody$expressed), ],
                 aes(x = wnnUMAP_1,
                     y = wnnUMAP_2),
                 colour = "gray", size = 0.1) +
      # cells with non-zero expression
      geom_point(data = colData_df_antibody[!is.na(colData_df_antibody$expressed), ],
                 aes(x = wnnUMAP_1,
                     y = wnnUMAP_2,
                     color = expressed,
                 ),size = 0.1)+
      scale_color_gradientn(
        name = "counts",
        na.value = "gray",
        colours = c("slateblue3", "royalblue1", "aquamarine3", "khaki", 383, "sienna1", "orangered4"),
        limits = c(1, max(3, upper_limit)),
        breaks = c(floor(upper_limit / 3), round(2 * (upper_limit / 3)), upper_limit)
      ) 
    print(antibody)
    print(paste0(my_ab," detected."))
    print(head(colData_df_antibody$expressed))
    antibody_plot <- ggplot(colData_df_antibody,
                            aes(x = wnnUMAP_1, y = wnnUMAP_2, color = expressed)) +
     # geom_point(data = colData_df_antibody[is.na(colData_df_antibody$expressed), ],
    #             size = rel(0.001),
    #             colour = "gray") +
    #  geom_point(data = colData_df_antibody[!is.na(colData_df_antibody$expressed), ],
      #geom_point(data = colData_df_antibody,
      geom_point(
                 #na.rm = TRUE,
                 size = 0.1) +
      #geom_point(data = colData_df_antibody[which(colData_df_antibody$expressed > max(thresholds[ab, ], upper_limit)), ],
      #           size = rel(0.001),
      #           colour = "orangered4") +
      scale_color_gradientn(
        name = "ln counts",
        na.value = "gray",
        colours = c("slateblue4", "royalblue1", "aquamarine3", "khaki", 383, "sienna1", "orangered4") ,
        limits = c(min_threshold, upper_limit),
        breaks = c(min_threshold,  medium_break, upper_limit),
        # breaks = c(min(thresholds[lookup[gene, ]$Antibody, ]), medium_break, upper_limit),
        labels = c(paste0("<", min_threshold),
                   medium_break,
                   upper_limit)
        #limits = c(thresholds[ab, ], upper_limit),
        #breaks = c(thresholds[ab, ], medium_break, upper_limit),
        #labels = c(
        #  paste("<", thresholds[ab, ], sep = ""),
        #  medium_break,
        #)
      ) +
      coord_fixed(ratio = 1) + theme_classic() +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(
          face = "bold", hjust = 0.5, vjust = 0.1, size = 8),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "pt"),
        legend.key.size = unit(0.8, "line"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.margin = margin(0.2, 0.2, 0.2, 0.2),
        axis.text = element_text(size = 5),
        aspect.ratio = 1
      ) +
      guides(fill = guide_legend(title = "counts"))
    print(antibody_plot)
  }
  antibody_plot
  # combine expression plot and antibody plot of current gene to one plot
  combined_plot  <- expr / antibody_plot
  plot_gene[[gene]] <- combined_plot
  ggplot2::ggsave(
    filename = paste("Test_RNAdataSlot", gene_name, "plot.png", sep = "_"),
    width = 21,
    height = 14,
    plot = plot_gene[[gene_name]])
  print("end")
}

# TODO What does this syntax mean?
layout <- c(
  area(1, 2, 5, 6),
  area(1, 7, 5, 11),
  area(1, 12, 5, 17),
#  area(1, 16, 5, 19),
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
numberOfPlots
# only plot 12 plots into one file
numberOfPages <- numberOfPlots %/% 12
if (numberOfPlots %% 12 != 0) {
  numberOfPages <-  numberOfPages + 1
}

# save plots
for (i in seq(1, numberOfPlots, 12)) {
  plot_collection <- p_final_ref + legend_ref  + p4# + cluster_umap + legend_cluster
  for (p in seq(0, 11)) {
    if (!(p+i>numberOfPlots)){
    plot_collection <- plot_collection + plot_gene[[p + i]]
    }
  }    
  plot_collection <- plot_collection + plot_layout(design = layout)
  final <- plot_collection #+ plot_annotation(caption = paste('Page ', (i + 12 - 1) / 12), " of ", numberOfPages) & theme(plot.caption = element_text(size = 10))
  out_expression_plots <- paste0(opt$outdir, opt$outName, "_")
  #out_expression_plots <- outName
  #out_expression_plots <- paste0("PBMC_TEST_")
  if (i + 11 > numberOfPlots) {
    ggplot2::ggsave(
      filename = paste(out_expression_plots, paste(lookup$Antibody[i ], "to", lookup$Antibody[numberOfPlots], "plot.png", sep = "_"), sep="."),
      width = 21,
      height = 14,
      plot = final
    )
  } else {
    ggplot2::ggsave(
      filename = paste(out_expression_plots, paste(lookup$Antibody[i ], "to", lookup$Antibody[i + 11], "plot.png", sep = "_"), sep="."),
      width = 21,
      height = 14,
      plot = final
    )
  }
}

# Dotplot

DefaultAssay(seurat_integrated) <- "adt"

my_adt <- t(seurat_integrated@assays$adt@counts)
meta <- as.data.frame(seurat_integrated@meta.data)
my_merge <- merge(my_adt,meta,by=0)
my_merge_to_compare <- my_merge 

head(my_merge)

unique(my_merge$sampleID)
for (ab in thresholds$Antibody){
  print(ab)
  for (sample  in unique(my_merge$sampleID)){
    my_merge[my_merge$sampleID==sample,paste0("ADT-",ab)] <- my_merge[my_merge$sampleID==sample,paste0("ADT-",ab)] - exp(thresholds[toupper(ab), sample])
    set_to_0 <- my_merge[my_merge$sampleID==sample,paste0("ADT-",ab)] <0
    my_merge[my_merge$sampleID==sample,paste0("ADT-",ab)][set_to_0] <-0
  }
  
}

# TODO
# Adapt this so it works for different panels & read in marker gene / protein list.
my_new_assay <- my_merge[,1:107]
dim(my_new_assay)
rownames(my_new_assay) <- my_new_assay$Row.names
my_new_assay$Row.names <- NULL
seurat_integrated@assays$adt_thresholds <-  CreateAssayObject(counts = t(my_new_assay))


Idents(seurat_integrated) <-seurat_integrated$celltype_final

DefaultAssay(seurat_integrated) <- "adt_thresholds"
DotPlot_adtThresholds <- DotPlot(seurat_integrated, features = c("ADT-CD45RA","ADT-CD45RO","ADT-CD45","ADT-CD19","ADT-CD20","ADT-CD22","ADT-CD73","ADT-CD1c","ADT-CD141",
                                                   "ADT-CD123","ADT-CD14","ADT-CD33","ADT-CD86","ADT-CD11b","ADT-KLRF1","ADT-CD16",
                                                   "ADT-CD57","ADT-CD3","ADT-CD4-RPA-T4","ADT-CD8a-RPA-T8","ADT-CD25"), 
                                 idents = c("B.cells","Dendritic.cells", "Dendritic.cells.Plasmacytoid","T.cells","T.cells.CD8","T.cells.CD4.naive", "NK.cells", "Monocytes"), 
                   assay="adt_thresholds", dot.scale = 10,scale = FALSE, scale.by = "size") + RotatedAxis()  + scale_colour_gradient2(low="gray", mid="lightblue", high="darkblue",na.value="darkblue")
  #guides(color = guide_colorbar(title = 'Scaled Average Expression')) 
DotPlot_adtThresholds
# TODO Adapt save
ggsave(filename = paste0(opt$outdir, opt$outName, "_DotPlot_ADT_Figure2_PBMC.png"),
       width = 40, height = 20, dpi = 600, units = "cm")


DotPlot_gex <- DotPlot(seurat_integrated, features = c("PTPRC","MS4A1","CD19","CD22","NT5E","CD1C","THBD","IL3RA","CD14","CD33","CD86","ITGAM","KLRF1","FCGR3A","B3GAT1","CD3E","CD4","CD8A","FOXP3","IL2RA"), 
                                 idents = c("B.cells","Dendritic.cells", "Dendritic.cells.Plasmacytoid","T.cells","T.cells.CD8","T.cells.CD4.naive", "NK.cells"), 
                                 assay="SCT", dot.scale = 10, col.min = 0,scale.by = "size") + RotatedAxis()  + scale_colour_gradient2(low="gray", mid="lightblue", high="darkblue",na.value="darkblue")+ 
  guides(color = guide_colorbar(title = 'Scaled Average Expression')) 
DotPlot_gex
  # TODO Adapt save
ggsave(filename = paste0(opt$outdir, opt$outName, "_DotPlot_GEX_Figure2_PBMC.png"),
       width = 40, height = 20, dpi = 600, units = "cm")

