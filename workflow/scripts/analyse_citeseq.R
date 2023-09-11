# Script to Analyse ADT Data in combination with GEXdata.
# Input: ADT Analysis on pooled Sample, GEX Analysis on individual Samples, lookup table matching gene and protein names, threshold table showing log thresholds per antibody
# Linda Grob


library(cowplot)
library(Seurat)
library(uwot)
library(SingleCellExperiment)
library(rjson)
library(stringr)
library(pheatmap)
library(patchwork)
library(ggplot2)
library(reshape)
library(ggridges)
library(dplyr)
library(tibble)
library(limma)
library(grid)
library(optparse)
library(rhdf5)
library(BREMSC)

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

option_list <- list(
  make_option("--RDS", type = "character", help = "Path to the atypical removed RDS file from the GEX pipline."),
  make_option("--cellrangerADT", type = "character", help = "Path to the cellranger ADT analysis folder (folder that contains feature, barcode and matrix file)."),
  make_option("--h5", type = "character", help = "Path to the hdf5 file of the most variable genes from the GEX pipline."),
  make_option("--colorConfig", type = "character", help = "Path to the colorConfig file used in the GEX pipeline."),
  make_option("--lookup", type = "character", help = "Path to a lookup table for gene & protein names."),
  make_option("--threshold", type = "character", help = "File containing individual log thresholds for each antibody."),
  make_option("--threads", type = "integer", help = "Number of threads that are available for the script. Recomended: 3-5"),
  make_option("--sampleName", type = "character", help = "SampleName. Needs to be exactly as used in the thresholds table."),
  make_option("--number_variable_genes", type = "integer", default = 500, help = "Number of variable genes that are included when calculating the UMAP embedding with RNA and ADT data."),
  make_option("--number_pca_adt", type = "integer", default = 20, help = "Number of PCA dimensions used when calculating the ADT UMAP. Cannot be larger than number of ADTs in experiment."),
  make_option("--output", type = "character", help = "Output directory.")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# set general options/parameters
options(stringsAsFactors = FALSE)

outfolder <- opt$output
dir.create(outfolder, showWarnings = FALSE)

outprefix <- paste(opt$output, opt$sampleName, sep = "/")
print("Files get the following prefix: ")
print(outprefix)

####
# Read in files
print("###################")
print("Reading in Input Files.")
print("###################")
print("Reading in RDS File.")
my_sce <- readRDS(opt$RDS)

print("Reading in cellranger ADT File.")
CellRangerADTfolder <- opt$cellrangerADT
print(CellRangerADTfolder)
CellRangerADT <- Read10X(data.dir = CellRangerADTfolder)

print("Read in h5 file containing variable genes used for umap clustering.")
vari <- h5read(opt$h5, "gene_attrs/gene_names")

print("Read in threshold file.")
ThresholdFile <- read.table(opt$threshold, header = TRUE, sep = "\t", row.names = 1)
thresholds <- ThresholdFile[opt$sampleName]
print("The following thresholds have been read for this sample:")
print(thresholds)

print(
  "Script assumes input thresholds to be logged. Raw Thresholds are then calculated to be the following:"
)
thresholds_exp <- round(exp(thresholds))
print(thresholds_exp)

print("Read in Color-Config")
config <- read.csv(opt$colorConfig, sep = "\t", stringsAsFactors = FALSE)
config$cell_type <- gsub(pattern = "([^_]+)_.*", "\\1", config$cell_type)
stopifnot(length(config$colour) == length(config$cell_type))
stopifnot(config$colour == unique(config$colour))
ct.color <- c(config$colour, "grey50", "black")
names(ct.color) <- c(config$cell_type, "uncertain", "unknown")
print(ct.color)
print("Read in Lookup Table")
lookup <- read.csv(opt$lookup, header = TRUE, sep = "\t")
print("The following lookup table is used to map gene names to antibodies and the other way round.")
print(lookup)

print("###################")
print("Start preprocessing:")
print("###################")
## This part is Sample Specific for the Levesque Pilot Samples. Should not be necessary later on!
print("Beautify rownames:")
rownames(x = CellRangerADT) <- gsub(
  pattern = "_TotalSeqC", replacement = "",
  x = rownames(x = CellRangerADT)
)
# Filter for identical cell barcodes in ADT & GEX analysis
print("Start filtering:")
print("Number of cells in GEX analysis before filtering:")
print(length(colnames(my_sce)))
print("Number of cells in ADT analysis before filtering:")
print(length(colnames(CellRangerADT)))
joint.bcs <- intersect(colnames(my_sce), colnames(CellRangerADT))
my_sce <- my_sce[, joint.bcs]
CellRangerADT <- CellRangerADT[, joint.bcs]
print("Number of cells in GEX analysis after filtering:")
print(length(colnames(my_sce)))
print("Number of cells in ADT analysis after filtering:")
print(length(colnames(CellRangerADT)))
# Check if a threshold can be assigned to every antibody.
print("Checking if a the antibodies are correctly named in the cellranger feature, the threshold and the lookup file.")
print(unique(rownames(CellRangerADT)))
print(unique(rownames(thresholds)))
stopifnot(setequal(unique(lookup$Antibody), unique(rownames(thresholds))))
print(unique(rownames(CellRangerADT)))
stopifnot(setequal(unique(lookup$Antibody), unique(rownames(CellRangerADT))))
# Defining list of antibodies and barcodes for future use in the script
print("Defining useful lists for future use:")
antibodies <- rownames(CellRangerADT)
antibodies <- antibodies[str_order(antibodies, numeric = TRUE)]
print("Alpha-numerical ordered list of antibodies used in ADT analysis.")
print(antibodies)
numberOfAntibodies <- length(antibodies)
print(paste("Found ", numberOfAntibodies, " Antibodies.", sep = ""))
barcodes <- colnames(CellRangerADT)
print(rownames(CellRangerADT))

# Define residual matrix on filtered SCE
residuals_matrix <- assay(my_sce, "pearson_resid")
residuals_matrix_variable <- residuals_matrix[rownames(residuals_matrix) %in% vari, ]

print("###################")
print("Creating summary tables:")
print("###################")
###
# Create a table containing phenograph_clusters, antibodies and celltypes per cellbarcode.
print("Create a table containing phenograph_clusters, antibodies and celltypes per cellbarcode.")
antibodies <- rownames(CellRangerADT)
barcodes <- colnames(CellRangerADT)
df_colData_cells <- as.data.frame(colData(my_sce), stringsAsFactors = FALSE)
table1 <- select(df_colData_cells, phenograph_clusters, celltype_final)
table1$antibody <- ""
for (i in seq_len(length(barcodes))) {
  table1$antibody[i] <- paste(antibodies[CellRangerADT[, i]>0], collapse = ",")
}
write.table(data.frame("barcodes" = rownames(table1), table1), file = paste(outprefix, "antibodies_and_celltypes_per_cell.rawcounts.tsv", sep = "."), sep = "\t", quote = FALSE)
print(paste("Saving it to ", outprefix, ".antibodies_and_celltypes_per_cell.rawcounts.tsv", sep = ""))

# Create a table containing fraction of cells associated with a Celltype expressing an antibody.
print("Create a table containing fraction of cells associated with a Celltype expressing an antibody.")

table2 <- as.data.frame(table(table1$celltype_final))
colnames(table2) <- c("Celltype", "NumberOfCells")
# Don't report celltypes with 0 cells.
table2 <- table2[which(table2$NumberOfCells > 0), ]
# Initiate empty antibody columns:
table2[, antibodies] <- 0

# Loop through celltypes & Antibodies
for (ab in antibodies) {
  for (ct in table2$Celltype) {
    occurence <- 0
    celltype_specific_abs <- table1[which(table1$celltype_final == ct), ]$antibody
    for (item in celltype_specific_abs) {
      if (grepl(paste(ab, ",", sep = ""), item) | grepl(paste(ab, "$", sep = ""), item)) {
        occurence <- occurence + 1
      }
    }
    fraction <- occurence / length(celltype_specific_abs)
    table2[which(table2$Celltype == ct), ab] <- round(fraction, 3)
  }
}
write.table(table2, file = paste(outprefix, "antibody_fraction_per_celltype.tsv", sep = "."), sep = "\t", quote = FALSE, row.names = FALSE)

print(paste("Saving it to ", outprefix, ".antibody_fraction_per_celltype.tsv", sep = ""))

# Prepare data for the fraction table

adtinfo <- as.data.frame(t(CellRangerADT))
adtinfo$barcodes <- colnames(CellRangerADT)
cell_attributes_extended <- plyr::join(as.data.frame(colData(my_sce)), adtinfo)
cell_attributes_extended <- cell_attributes_extended[complete.cases(cell_attributes_extended), ]

# Create a table containing fraction of cells expressing the corresponding gene to an antibody.
print("Create a table containing fraction of cells expressing the corresponding gene to an antibody.")
count_matrix <- assay(my_sce, "normcounts")
table3 <- lookup
table3$number_cells_per_AB <- 0
table3$fraction_total_cells_with_ab <- 0
table3$genes_expressed <- 0
table3$genes_not_expressed <- 0
table3$fraction_ab_labeled_cells_expressing_gene <- 0
colnames(table3) <- c("Gene", "Antibody", "number_cells_per_AB", "fraction_total_cells_with_ab", "genes_expressed", "genes_not_expressed", "fraction_ab_labeled_cells_expressing_gene")
thresholds[is.na(thresholds)] <- 0
for (ab in seq(length(lookup$Antibody))) {
  count_expressed <- 0
  count_not_expressed <- 0
  count_ab <- sum(cell_attributes_extended[, lookup[ab, ]$Antibody]>thresholds[lookup[ab, ]$Antibody, ])
  for (index in seq(length(cell_attributes_extended$barcodes))) {
    barcode <- cell_attributes_extended$barcodes[index]
    gene <- lookup[ab, ]$Gene
    if (is.na(match(gene, rownames(count_matrix)))) {
      next
    }
    if (cell_attributes_extended[index, lookup[ab, ]$Antibody] > thresholds[lookup[ab, ]$Antibody, ]){
      if (count_matrix[gene, barcode] != 0) {
        count_expressed <- count_expressed + 1
      }
      else if (count_matrix[gene, barcode] == 0 ) {
        count_not_expressed <- count_not_expressed + 1
      }
    }
  }
  table3[ab, ]$genes_not_expressed <- count_not_expressed
  table3[ab, ]$genes_expressed <- count_expressed
  table3[ab, ]$fraction_ab_labeled_cells_expressing_gene <- round(count_expressed / (count_ab), 2)
  table3[ab, ]$number_cells_per_AB <- count_ab
  table3[ab, ]$fraction_total_cells_with_ab <- round(count_ab / length(cell_attributes_extended$barcodes), 2)
}
print(paste("Saving it to ", outprefix, ".fraction_of_genes_expressed_per_antibody.tsv", sep = ""))
write.table(table3, file = paste(outprefix, "fraction_of_genes_expressed_per_antibody.tsv", sep = "."), sep = "\t", row.names = FALSE, quote = FALSE)

print("###################")
print("Start plotting:")
print("###################")
###   Plot antibody heatmap   ###
matrix_CellRangerADT <- as.matrix(CellRangerADT)
log_matrix_CellRangerADT <- log(as.matrix(CellRangerADT) + 1)
matrices_list <-
  list(log_matrix_CellRangerADT, matrix_CellRangerADT)
names(matrices_list) <- c("log", "raw")
# get colours for heatmap
col_palette <- viridisLite::plasma(255)

# get data for cell type bar in heatmap
final_celltypes <- as.data.frame(colData(my_sce)$celltype_final)
# increase length of legend title to make sure the whole legend is displayed
additional_ws <- max(nchar(names(ct.color)))
celltype_string <- paste0("Cell type", strrep(" ", additional_ws))
names(final_celltypes) <- celltype_string
rownames(final_celltypes) <- colData(my_sce)$barcodes
stopifnot(rownames(final_celltypes) == colnames(matrix_CellRangerADT))

# use colours from colour config for cell type bar
list_colours <- list(ct.color)
names(list_colours) <- celltype_string

# plot two heatmaps once using log, and once using raw counts
for (hm_version in c("log", "raw")) {
  print(paste0("Plot heatmap of ", hm_version, " ADT counts"))
  # get file name for heatmap
  filename_hm_phmp <-
    paste0(outprefix, ".", hm_version, "_antibody_counts_heatmap.png")
  print(filename_hm_phmp)
  hm_phmp <- pheatmap(
    matrices_list[[hm_version]],
    scale = "none",
    clustering_method = "ward.D2",
    show_colnames = FALSE,
    color = col_palette,
    cluster_rows = FALSE,
    annotation_names_col = T,
    fontsize = 11,
    main = paste0("Antibody counts ", hm_version),
    treeheight_col = 45,
    annotation_col = final_celltypes,
    annotation_colors = list_colours
  )
  # save heatmap
  ggsave(
    filename_hm_phmp,
    hm_phmp,
    width = 40,
    height = max(20, round(length(antibodies) * 0.4)),
    units = "cm",
    dpi = 300
  )
}

###
# Ridgeplot ABs
###
print("Plotting Ridgeplot Antibodies:")
ridgeout <- paste(outfolder, "/RidgeplotsSingleSample/", sep = "")
dir.create(ridgeout)

plot_list <- list()
for (ab in antibodies) {
  print(paste("Working on ", ab, " using threshold ", thresholds[ab, ], ".", sep =
                ""))
  abcounts <- as.matrix(CellRangerADT[ab, ])
  colnames(abcounts) <- ab
  abcounts <- melt(abcounts)
  colnames(abcounts) <- c("barcode", "antibody", "value")
  # To avoid issues with zero counts add +1 to all antibody counts before talking the log
  abcounts$value <- abcounts$value + 1
  abcounts$ln_umi <- log(abcounts$value)
  # Treat special case, if no count has been detected.
  if (all(is.na(abcounts$ln_umi))) {
    abPlot <- ggplot() +
      theme_void() +
      theme(plot.title = element_text(vjust = -50, hjust = 20)) +
      ggtitle(paste(ab, ": not detected"))
  } else {
    abPlot <-
      ggplot(abcounts,
             aes(
               x = ln_umi,
               y = antibody,
               fill = antibody,
               height = ..ndensity..
             )) +
      geom_density_ridges() +
      scale_x_continuous(limits = c(0, max(abcounts$ln_umi) + 1)) +
      theme_ridges() +
      theme(
        legend.position = "none",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm")
      )
    # Draw a line showing threshold if threshold is > 0
    if (as.numeric(thresholds[ab, ]) > 0) {
      abPlot <-
        abPlot + geom_vline(xintercept = as.numeric(thresholds[ab, ]))
    }
  }
  plot_list[[ab]] <- abPlot
}
print("Saving ridgeplots.")
numberOfPlots <- length(plot_list)
for (i in seq(1, numberOfPlots, 12)) {
  plot_sub <- plot_list[i:min(i + 11, numberOfPlots)]
  plot_group <-
    plot_grid(plotlist = plot_sub,
              ncol = 4,
              align = "v")
  if (i + 11 > numberOfPlots) {
    ggplot2::ggsave(
      filename = paste(
        paste(ridgeout, opt$sampleName, ".Antibody", sep = ""),
        antibodies[i],
        "to",
        antibodies[numberOfPlots],
        "ridgeplot.png",
        sep = "_"
      ),
      width = 21,
      height = 14,
      plot = plot_group
    )
  } else {
    ggplot2::ggsave(
      filename = paste(
        paste(ridgeout, opt$sampleName, ".Antibody", sep = ""),
        antibodies[i],
        "to",
        antibodies[i + 11],
        "ridgeplot.png",
        sep = "_"
      ),
      width = 21,
      height = 14,
      plot = plot_group
    )
  }
}
###
# Plot Antibodies per Celltype
# Since this are potentially many plots, create an extra subfolder for the output.
print("Plotting Antibodies per Celltype and Calculating table:")
abcelltypeout <-
  paste(outfolder, "/antibodies_per_celltype_plots/", sep = "")
dir.create(paste(abcelltypeout), showWarnings = FALSE)

# Create a dataframe containing adt and gex information per barcode
df_colData_cells <- as.data.frame(colData(my_sce), stringsAsFactors = FALSE)
CellRangerADTextended <- rownames_to_column(as.data.frame(t(CellRangerADT)))
colnames(CellRangerADTextended)[1] <- "barcodes"
cell_attributes <- plyr::join(df_colData_cells, CellRangerADTextended, by = "barcodes")
print("Saving cell data table, containing the following columns:")
print(colnames(cell_attributes))
write.table(cell_attributes, paste(outprefix, ".cell_attributes.tsv", sep = ""),
            quote = FALSE, row.names = FALSE, sep = "\t")


# Save SCE object with GEX results AND cellranger ADT results in colData table
cat("\n\nSaving RDS file containing the whole GEX SCE object and the cellranger ADT columns.\n")
stopifnot(CellRangerADTextended$barcodes == colData(my_sce)$barcodes)
for (column in 1:length(CellRangerADTextended)) {
  columnName <- names(CellRangerADTextended)[column]
  columnName <- paste0("ADT_", columnName)
  cat("Add column", columnName, "to SCE object\n")
  colData(my_sce)[columnName] <- CellRangerADTextended[[column]]
}
saveRDS(my_sce, paste(outprefix, ".GEX_cellrangerADT_SCE.RDS", sep = ""))


id.final.ct <- match(levels(as.factor(cell_attributes$celltype_final)), names(ct.color))
id.major.ct <- match(levels(as.factor(cell_attributes$celltype_major)), names(ct.color))
celltype_list <- list(id.final.ct, id.major.ct)
names(celltype_list) <- c("celltype_final", "celltype_major")

# Save tables showing proportions of celltypes.
print("Saving tables containing proportions of celltype (major)")
write.table(
  round(prop.table(table(
    cell_attributes$celltype_major
  )), 2),
  file = paste0(outprefix, ".major_celltype.fractions.txt"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

# plot ridgeplots for each ADT, per major/minor cell types
plot_list <- list()
for (type in c("celltype_final", "celltype_major")) {
  print(paste("Working on plotting ridgeplots using ", type, ".", sep = ""))
  for (ab in antibodies) {
    print(paste("Plotting ", ab, sep = ""))
    plot_data <- select(cell_attributes, all_of(type))
    plot_data$antibody <- cell_attributes[, ab]
    colnames(plot_data) <- c(type, "antibody")
    plot_data$antibody <- plot_data$antibody + 1
    plot_data$ln_umi <- log(plot_data$antibody)
    plot_data <- plot_data[complete.cases(plot_data), ]
    s <- ggplot(plot_data,
             aes_string(
               x = "ln_umi",
               y = type,
               fill = type,
               height = "..ndensity.."
             )) +
      geom_density_ridges(
        alpha = 0.82,
        jittered_points = TRUE,
        point_alpha = 0.4,
        point_size = 0.4
      ) +
      theme_ridges() +
      scale_x_continuous(limits = c(0, max(plot_data$ln_umi) + 1)) +
      ggtitle(ab) +
      guides(fill = FALSE) +
      scale_fill_manual(name = "Cell type",
                        values = ct.color[unlist(celltype_list[type])],
                        drop = F) +
      theme(axis.title.y = element_blank(),
            plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"))
    # Add line indicating threshold if threshold is larger than 0
    if (as.numeric(thresholds[ab, ]) > 0) {
      s <- s + geom_vline(xintercept = as.numeric(thresholds[ab, ]))
    }
    plot_list[[ab]] <- s
  }
  numberOfPlots <- length(plot_list)
  print("Saving celltype ridgeplots.")
  for (i in seq(1, numberOfPlots, 6)) {
    plot_sub <- plot_list[i:min(i + 5, numberOfPlots)]
    plot_group <-
      plot_grid(plotlist = plot_sub,
                ncol = 3,
                align = "v")
    if (i + 5 > numberOfPlots) {
      ggplot2::ggsave(
        filename = paste(
          paste(abcelltypeout, opt$sampleName, ".Antibody.", type, sep = ""),
          antibodies[i],
          "to",
          antibodies[numberOfPlots],
          "ridgeplot.png",
          sep = "_"
        ),
        width = 20,
        height = 14,
        plot = plot_group
      )
    } else {
      ggplot2::ggsave(
        filename = paste(
          paste(abcelltypeout, opt$sampleName, ".Antibody.", type, sep = ""),
          antibodies[i],
          "to",
          antibodies[i + 5],
          "ridgeplot.png",
          sep = "_"
        ),
        width = 20,
        height = 14,
        plot = plot_group
      )
    }
  }
}


print("###################")
print("Start clustering:")
print("###################")
set.seed(184)
# Coulours used later on for clustering
cols33 <- c("red2", "green4", "blue2", "cyan2", "yellow1", "purple", "brown",
            "chocolate1", "chartreuse2", "darkgoldenrod3", "steelblue1", "slateblue3", "olivedrab4", "gold2",
            "violetred3", "darkcyan", "orchid3", "darksalmon", "darkslategrey", "khaki", "indianred2" ,"magenta", "slategray2",
            "olivedrab1", "mediumaquamarine", "hotpink", "yellow3",
            "bisque4", "darkseagreen1", "dodgerblue3",
            "deeppink4", "sienna4", "mediumorchid4")
print("Defined colours for clustering; only 33 available, if more clusters are found, this might crash.")
### Cluster
print("Start clustering on ADT and GEX data.")
RNA_count_data <- assay(my_sce, "counts")
# reduce RNA count data to number of most variable genes specified as command line argument
subset_vari_rowData <- rowData(my_sce)[rownames(rowData(my_sce)) %in% vari, ]
subset_vari_ordered <- subset_vari_rowData[order(subset_vari_rowData$residual_variance, decreasing = TRUE), ]
vari_ordered <- subset_vari_ordered$gene_names
RNA_count_data <- RNA_count_data[rownames(RNA_count_data) %in% vari_ordered[1:opt$number_variable_genes], ]
# Clustering taking the number of clusters from the phenograph)
cluster <- BREMSC(CellRangerADT, RNA_count_data, K = max(cell_attributes$phenograph_clusters), nChains = opt$threads, nMCMC = 1000)
cluster.df <- as.data.frame(cluster$clusterID)
cluster.df$barcodes <- colnames(CellRangerADT)
resort_clusterIDs <- as.data.frame(sort(table(cluster$clusterID), decreasing = TRUE))
colnames(resort_clusterIDs) <- c("oldID", "freq")
resort_clusterIDs$newIDs <- rownames(resort_clusterIDs)
cluster.df$newIDs <- resort_clusterIDs$newIDs[match(unlist(cluster$clusterID), resort_clusterIDs$oldID)]

###
# Prepare data for the main plot showing gex vs. ab expression and the third table
# get UMAP coordinates into data frame
dir.create(paste(outfolder, "/ExpressionPlots/", sep = ""))
dir.create(paste(outfolder, "/ExpressionPlots/gex_based_embedding/", sep = ""))
dir.create(paste(outfolder, "/ExpressionPlots/adt_based_embedding/", sep = ""))

print("###################")
print("Start creating Expression Plots:")
print("###################")

for (type in c("adt", "gex", "adt_gex")) {
  if (type == "adt") {
    cell_attributes <- plyr::join(df_colData_cells, CellRangerADTextended, by = "barcodes")
    print("Working on Adt based UMAP embedding.")
    combout <- paste(outfolder, "/ExpressionPlots/adt_based_embedding/", opt$sampleName, ".ADT_only", sep = "")
    umap.adt <- umap(t(as.matrix(CellRangerADT)), n_neighbors = 30, pca = opt$number_pca_adt, spread = 1, min_dist = 0.3, ret_nn = T)
    reducedDim(my_sce, "umap_adt") <- umap.adt$embedding
    metadata(my_sce) <- c(metadata(my_sce), list(umap_adt = umap.adt$nn$euclidean))
    umap_coord <- as.data.frame(reducedDims(my_sce)$umap_adt)
    names(umap_coord) <- c("umap1", "umap2")
    umap_coord$barcodes <- barcodes
    # make sure all cell dimensions are the same
    stopifnot(my_sce$barcodes == colData(my_sce)$barcodes)
    stopifnot(my_sce$barcodes == umap_coord$barcodes)
    # add umap coordinates to the meta data table of cells
    cell_attributes <- plyr::join(cell_attributes, umap_coord)
  } else if (type == "gex") {
    cell_attributes <- plyr::join(df_colData_cells, CellRangerADTextended, by = "barcodes")
    print("Working on Gex based UMAP embedding.")
    combout <- paste(outfolder, "/ExpressionPlots/gex_based_embedding/", opt$sampleName, ".GEX_only", sep = "")
    # get UMAP coordinates into data frame
    umap_coord <- as.data.frame(reducedDims(my_sce)$umap_hvg)
    names(umap_coord) <- c("umap1", "umap2")
    umap_coord$barcodes <- rownames(umap_coord)
    # make sure all cell dimensions are the same
    stopifnot(my_sce$barcodes == colData(my_sce)$barcodes)
    stopifnot(my_sce$barcodes == umap_coord$barcodes)
    # add umap coordinates to the meta data table of cells
    cell_attributes <- plyr::join(cell_attributes, umap_coord)
  } else {
    cell_attributes <- plyr::join(df_colData_cells, CellRangerADTextended, by = "barcodes")
    combout <- paste(outfolder, "/ExpressionPlots/", opt$sampleName, ".Combined_ADT_GEX", sep = "")
    print("Working on joined Adt and Gex based UMAP embedding.")
    CellRangerreduced  <-  CellRangerADT[, colnames(residuals_matrix_variable)] + 1
    residuals_matrix_variable_reduced <- residuals_matrix_variable[, colnames(CellRangerADT)]
    CombinedADTGEX <- rbind(residuals_matrix_variable_reduced, CellRangerreduced)
    umap.adt <- umap(t(CombinedADTGEX), scale = TRUE, n_neighbors = 30, pca = 50, spread = 1, min_dist = 0.4, ret_nn = T)
    # get UMAP coordinates into data frame
    reducedDim(my_sce, "umap_adt") <- umap.adt$embedding
    metadata(my_sce) <- c(metadata(my_sce), list(umap_adt = umap.adt$nn$euclidean))
    #umap_coord <- as.data.frame(reducedDims(my_sce)$umap_hvg)
    umap_coord <- as.data.frame(reducedDims(my_sce)$umap_adt)
    names(umap_coord) <- c("umap1", "umap2")
    umap_coord$barcodes <- rownames(umap_coord)
    # add umap coordinates to the meta data table of cells
    cell_attributes <- plyr::join(cell_attributes, umap_coord)
  }

  # Plot Clustering
  cluster_attributes <- plyr::join(cell_attributes, cluster.df)
  cluster_attributes$cluster <- as.integer(cluster$clusterID)
  cluster_attributes$newIDs <- as.integer(cluster_attributes$newIDs)

  cluster_umap <- ggplot(cluster_attributes, aes(x = umap1, y = umap2, color = as.factor(newIDs))) +
    geom_point(size = 1) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 15)) +
    scale_color_manual(name = "BREMSC", values = cols33)+
    theme(aspect.ratio = 1) +
    theme(legend.position = "none")

  cluster_umap_final <- ggplot(cluster_attributes, aes(x = umap1, y = umap2, color = as.factor(newIDs))) +
    geom_point(size = 1) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 15)) +
    scale_color_manual(name = "BREMSC", values = cols33)+
    theme(aspect.ratio = 1) +
    theme(
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 10),
      legend.key.width = unit(0.2, "line"),
      legend.key.height = unit(0.01, "line"),
      legend.spacing.y = unit(0.1, "line"),
      legend.spacing.x = unit(0.1, "line")
  )

  legend_cluster <- cowplot::get_legend(cluster_umap_final)

  print("Antibody / Gene Expression - Plots are created.")
  gene.list.all <- unique(lookup$Gene)
  gene.list.all <- data.frame(Gene = unique(lookup$Gene))

  id.first.ct <- match(levels(cell_attributes$celltype_major), names(ct.color))
  id.final.ct <- match(levels(cell_attributes$celltype_final), names(ct.color))

  # have matrix with normcounts
  # remove all genes that have 0 counts in all cells
  normcounts_all.zero.removed <- normcounts(my_sce)
  mask_all_zero <- apply(normcounts_all.zero.removed, 1, sum) > 0
  normcounts_all.zero.removed <-
  normcounts_all.zero.removed[mask_all_zero, ]
  stopifnot(length(normcounts_all.zero.removed[, 1]) == sum(apply(normcounts_all.zero.removed, 1, sum) > 0))

  # keep only genes that are in matrix
  gene.list <- gene.list.all[gene.list.all$Gene %in% rownames(normcounts_all.zero.removed), ]

  # sort gene names alphabetically for better overview in plot
  gene.list <- sort(gene.list)
  gene.list.all <- sort(gene.list.all$Gene)


  # get indices of genes in matrix
  list_groups <- lapply(seq_along(gene.list), function(x) {
    indices <- match(gene.list[[x]], rownames(normcounts_all.zero.removed))
    indices
  })

  # Define some general used graphical parameters
  fontsize <- theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10))
  theme_set(theme_bw(4) + fontsize)

  # plot UMAP (based on highly variable genes), colours = final cell type, as reference for expression plots
  p_final_ref <-
    ggplot(cell_attributes, aes(x = umap1, y = umap2, color = celltype_final)) +
    scale_color_manual(name = "Cell type",
                     values = ct.color[id.final.ct],
                     drop = F) +
    geom_point(size = 1) +
    theme(aspect.ratio = 1) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    theme(legend.position = "none")

  p_legend_ref <-
    ggplot(cell_attributes, aes(x = umap1, y = umap2, color = celltype_final)) +
    scale_color_manual(name = "Cell type",
                     values = ct.color[id.final.ct],
                     drop = F) +
    geom_point(size = 1) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 15)) +
    theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 10),
    legend.key.width = unit(0.2, "line"),
    legend.key.height = unit(0.01, "line"),
    legend.spacing.y = unit(0.1, "line"),
    legend.spacing.x = unit(0.1, "line")
    )

  legend_ref <- cowplot::get_legend(p_legend_ref)

  # plot UMAP (based on highly variable genes), colours = final cell type, as reference for expression plot
  # get all genes that have zero counts in all cells
  mask.not.expressed <-
    !(gene.list.all %in% rownames(normcounts_all.zero.removed[unlist(list_groups), , drop = F]))
  genes.not.expressed <- sort(gene.list.all[mask.not.expressed])

  # get matrix with all cells and all genes of the respective group
  matrix_group <-
  normcounts_all.zero.removed[unlist(list_groups), ,  drop = F]

  # check first if any of the group's genes are expressed in this sample
  lookup <- lookup[str_order(lookup$Antibody, numeric = TRUE), ]
  lookup_expressed <- lookup[-which(lookup$Gene %in% genes.not.expressed), ]
  lookup_not_expressed <- lookup[which(lookup$Gene %in% genes.not.expressed), ]
  rownames(lookup_expressed) <- c(seq_len(length(lookup_expressed$Gene)))
  rownames(lookup_not_expressed) <- c(seq_len(length(lookup_not_expressed$Gene)))

  # get empty list that will be filled with expression plots per gene
  plot_gene <- list()

  expout <- paste(outfolder, "/ExpressionPlots/", opt$sampleName, sep = "")
  # loop over the expressed genes of the current group
  for (gene in seq(length(lookup$Gene))) {
    merged_plots <- list()
    gene_name <- lookup$Gene[gene]
    print(paste("Working on ", gene_name, ".", sep = ""))
    # combine cell_attributes with expression of current gene
    cell_attributes_gene <- cell_attributes
    if (gene_name %in% rownames(matrix_group)) {
      cell_attributes_gene$normcounts <- matrix_group[gene_name, ]
      # maximum count found for current gene
      max_count <- max(cell_attributes_gene$normcounts)
      # upper limit of gene expression colour scale is either maximum count or 3
      upper_limit <- max(3, max_count)
      # Plotting only gene and antibody name if they are not the same
      if (toupper(gene_name) == toupper(lookup[gene, ]$Antibody)) {
        legend <- gene_name
      } else {
        legend <- paste(gene_name, " / ", lookup[gene, ]$Antibody, sep = "")
      }
      # plot gene expression
      expr <-
        ggplot(cell_attributes_gene, aes(x = umap1, y = umap2)) +
        geom_point(colour = "gray", size = rel(0.001)) +
        geom_point(data = cell_attributes_gene[which(cell_attributes_gene$normcounts == 0), ],
                 size = rel(0.001),
                 colour = "slateblue4") +
        geom_point(data = cell_attributes_gene[which(cell_attributes_gene$normcounts > 0), ],
                 aes(color = normcounts),
                 size = rel(0.001)) +
        scale_color_gradientn(
          name = "counts",
          na.value = "gray",
          colours = c(
          "slateblue3",
          "royalblue1",
          "aquamarine3",
          "khaki",
          383,
          "sienna1",
          "orangered4"
        ),
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
          face = "bold",
          hjust = 0.5,
          vjust = 0.1,
          size = 14
        ),
        plot.margin = margin(1, 1, 1, 1, "pt"),
        legend.key.size = unit(0.8, "line"),
        legend.text = element_text(size = rel(3)),
        legend.title = element_text(size = rel(3)),
        legend.margin = margin(0.2, 0.2, 0.2, 0.2),
        axis.text = element_text(size = rel(3)),
        aspect.ratio = 1
      )
    } else {
      cell_attributes_gene$normcounts <- matrix_group[1, ]
      cell_attributes_gene$normcounts <- 0
      # maximum count found for current gene
      max_count <- max(cell_attributes_gene$normcounts)
      # upper limit of gene expression colour scale is either maximum count or 3
      upper_limit <- max(3, max_count)
      if (toupper(gene_name) == toupper(lookup[gene, ]$Antibody)) {
        legend <- gene_name
      } else {
        legend <- paste(gene_name, " / ", lookup[gene, ]$Antibody, sep = "")
      }
      # plot gene expression
      expr <-
        ggplot(cell_attributes_gene, aes(x = umap1, y = umap2)) +
        geom_point(colour = "gray", size = rel(0.001)) +
        geom_point(data = cell_attributes_gene[which(cell_attributes_gene$normcounts == 0), ],
                   size = rel(0.001),
                   colour = "gray") +
        scale_color_gradientn(
          name = "counts",
          na.value = "gray",
          colours = c(
            "slateblue3",
            "royalblue1",
            "aquamarine3",
            "khaki",
            383,
            "sienna1",
            "orangered4"
          ),
          limits = c(1, max(3, upper_limit)),
          breaks = c(floor(upper_limit / 3), round(2 * (upper_limit / 3)), upper_limit)
        ) +
        coord_fixed(ratio = 1) +
        ggtitle(legend) +
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
            size = 14
          ),
          plot.margin = margin(1, 1, 1, 1, "pt"),
          legend.key.size = unit(0.8, "line"),
          legend.text = element_text(size = rel(3)),
          legend.title = element_text(size = rel(3)),
          legend.margin = margin(0.2, 0.2, 0.2, 0.2),
          axis.text = element_text(size = rel(3)),
          aspect.ratio = 1
        )
    }
    # combine cell_attributes with antibody
    cell_attributes_antibody <- cell_attributes
    cell_attributes_antibody$expressed <- as.numeric(cell_attributes[, lookup[gene, ]$Antibody]) + 1
    cell_attributes_antibody$expressed <- round(log(cell_attributes_antibody$expressed))
    ab <- lookup[gene, ]$Antibody

        # plot antibody plot
    upper_limit <- max(thresholds[ab, ] * 3, max(cell_attributes_antibody$expressed), 3)
    medium_break <- (thresholds[ab, ] + upper_limit) / 2
    upper_break <- (medium_break + upper_limit) / 2
    lower_break <- (medium_break + thresholds[ab, ]) / 2
    if (all(cell_attributes_antibody$expressed == 0)) {
      antibody_plot <- ggplot(cell_attributes_antibody,
        aes(x = umap1, y = umap2, color = expressed)) +
        geom_point(data = cell_attributes_antibody[which(cell_attributes_antibody$expressed == 0), ],
                 size = rel(0.001),
                 colour = "gray") +
        scale_color_gradientn(
          name = "ln counts",
          na.value = "slateblue4",
          colours = c(
            "slateblue4",
            "royalblue1",
            "aquamarine3",
            "khaki",
            383,
            "sienna1",
            "orangered4"
          ),
          limits = c(thresholds[ab, ], upper_limit),
          breaks = c(thresholds[ab, ], medium_break, upper_limit),
          labels = c(
            paste("<", thresholds[ab, ], sep = ""),
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
    } else {
      antibody_plot <- ggplot(cell_attributes_antibody,
                            aes(x = umap1, y = umap2, color = expressed)) +
        geom_point(data = cell_attributes_antibody[which(cell_attributes_antibody$expressed <= thresholds[ab, ]), ],
                 size = rel(0.001),
                 colour = "gray") +
        geom_point(data = cell_attributes_antibody[which(cell_attributes_antibody$expressed > thresholds[ab, ]), ],
                 na.rm = TRUE,
                 size = rel(0.001)) +
        geom_point(data = cell_attributes_antibody[which(cell_attributes_antibody$expressed > max(thresholds[ab, ], upper_limit)), ],
                 size = rel(0.001),
                 colour = "orangered4") +
        scale_color_gradientn(
          name = "ln counts",
          na.value = "slateblue4",
          colours = c(
            "gray",
            "slateblue4",
            "royalblue1",
            "aquamarine3",
            "khaki",
            383,
            "sienna1",
            "orangered4"
          ),
          limits = c(thresholds[ab, ], upper_limit),
          breaks = c(thresholds[ab, ], medium_break, upper_limit),
          labels = c(
            paste("<", thresholds[ab, ], sep = ""),
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
    }
    # combine expression plot and antibody plot of current gene to one plot
    plot_gene[[gene]] <- expr / antibody_plot
  }




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


  print("Saving Expression Plots:")
  numberOfPlots <- length(plot_gene)
  numberOfPages <- numberOfPlots %/% 12
  if (numberOfPlots %% 12 != 0) {
    numberOfPages <-  numberOfPages + 1
  }
  for (i in seq(1, numberOfPlots, 12)) {
    plot_collection  <- p_final_ref + legend_ref + cluster_umap + legend_cluster
    for (p in seq(1, 12)) {
      plot_collection <- plot_collection + plot_gene[p + i]
    }
    plot_collection <- plot_collection + plot_layout(design = layout)
    final <- plot_collection + plot_annotation(caption = paste('Page ', (i + 12 - 1) / 12), " of ", numberOfPages) & theme(plot.caption = element_text(size = 10))
    if (i + 12 > numberOfPlots) {
      ggplot2::ggsave(
        filename = paste(combout, lookup$Antibody[i + 1], "to", lookup$Antibody[numberOfPlots], "plot.png", sep = "_"),
        width = 21,
        height = 14,
        plot = final
      )
    } else {
      ggplot2::ggsave(
        filename = paste(combout, lookup$Antibody[i + 1], "to", lookup$Antibody[i + 12], "plot.png", sep = "_"),
        width = 21,
        height = 14,
        plot = final
      )
    }
  }

}
