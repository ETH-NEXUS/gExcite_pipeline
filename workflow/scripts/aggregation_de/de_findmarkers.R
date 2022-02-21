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
## optparse 1.7.1

## GENERAL:
## This script allows for Differential expression analysis to be directly performed on an integrated Seurat object, without the need to access single samples or clinical data files.
## Please refer to the readme "doc/readme_differential_expression" for details about the DE step, differences between the available methods and suggestions on to which mehod to choose depending on your data
## Several plots are produced alongside the actual output to provide insight on the data produced

### ADDITIONAL IMPORTANT INFORMATION:
## This script takes as input a Seurat object that must contain count data for a single cell experiment.
## Contrasts are computed based on a metadata variable of interest that must be available in the Seurat object.
## For a given feature selected for the contrasts, the script performs DE on all possible pairwise contrasts. In addition also the contrasts <value> vs all is performed
## If the tool gives a warning about the "limma" package missing, it is strongly encouraged to stop the script, install "limma" and restart the DE analysis. The presence of the limma package drastically cuts the run times
## All options are required and must have a value set. In the future there might be a migration away from optparse and to a package allowing to
## automatically check if a required option has a value or not.
## The heatmap will be generated using a different script

## INPUT:
## - The path where to find the RDS file storing the Seurat object containing the data
## - The feature used to compute contrasts. It must be an existing metadata column of the provided Seurat object. The script will calculate all possible pairwise contrasts based on the content of the feature
## - The feature used to split the Seurat object. It must be an existing metadata column of the provided Seurat object. If set to NULL, all Seurat object will be used, otherwise, the script will subset the seurat object based on the feature values and perform DE every single one separately
## - The assay slot of the Seurat object from which to retrieve the expression
## - The directory in which the output must be stored
## - The FDR cutoff to use to define a gene as significant
## - The fold-change cutoff to use to define a gene as strong

## OUTPUT:
## - For each element of the grouping feature (e.g. for each celltype if grouping is done by celltype)
##   AND for each contrast generated based on the test variable:
##   - Volcano plot
##   - RDS object containing the Table with the full results of the DE analysis
##   - Tab-delimited table containing the filtered results of the DE analysis

suppressPackageStartupMessages({
  library(ggplot2)
  library(Seurat)
  library(optparse)
})

option_list <- list(
  make_option("--seuratFile", type = "character", help = "Path to the RDS file containing the Seurat object with the counts and metadata"),
  make_option("--contrastFeature", type = "character", help = "Name of the feature to use to compute the contrasts. It must be an existing metadata column in the Seurat object"),
  make_option("--splitFeature", type = "character", default = "", help = "Name of the feature to use to split the Seurat object. It must be an existing metadata column in the Seurat object"),
  make_option("--assay", type = "character", help = "Assay slot of the Seurat object to use for the analysis"),
  make_option("--outdir", type = "character", "Path to output directory"),
  make_option("--fdr.cut", type = "double", default = 0.05, help = "The FDR cutoff for the DE analysis"),
  make_option("--fc.cut", type = "double", default = 1.5, help = "The fold-change cutoff for the DE analysis"),
  make_option("--test.use", type = "character", default = "wilcox", help = "Test to use for the DE analysis. Please use 'LR' or 'MAST' is you wish to add covariates/confounding factors"),
  make_option("--latent.vars", type = "character", default = NULL, help = "single string or array of covariates/confounding factors. Only meaningful if 'test.use' is set to either 'LR' or 'MAST'")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

"%&%" <- function(a, b) paste(a, b, sep = "")

seuratFile <- opt$seuratFile
contrastFeature <- opt$contrastFeature
splitFeature <- opt$splitFeature
myassay <- opt$assay
outdir <- opt$outdir
myfdr <- opt$fdr.cut
myfc <- opt$fc.cut
test.use <- opt$test.use
latent.vars <- opt$latent.vars

message("Loading the Seurat object.\nDepending on the size of project, this may take a while and use a lot of memory...")
seurat_integrated <- readRDS(seuratFile)


# Input sanity check
if (splitFeature != "") {
  # Gracefully die if the feature to split for is not in the Seurat object
  if (!any(colnames(seurat_integrated@meta.data) %in% splitFeature)) {
    stop("ERROR: The selected split feature is not in the Seurat object")
  }
}
# Gracefully die if the feature from which to compute contrasts is not in the Seurat object
if (!any(colnames(seurat_integrated@meta.data) %in% contrastFeature)) {
  stop("ERROR: The selected feature from which to compute contrasts is not in the Seurat object")
}
# or if it has only one value
if (length(unique(seurat_integrated@meta.data[, which(colnames(seurat_integrated@meta.data) == contrastFeature)])) == 1) {
  stop("ERROR: The feature from which to compute contrasts only contains one level. No contrast is possible")
}
# Gracefully die if the assay does not exist
if (!any(names(seurat_integrated@assays) %in% myassay)) {
  stop("ERROR: The selected assay is not in the Seurat object")
}

contrastFeaturePos <- which(colnames(seurat_integrated@meta.data) == contrastFeature)
splitFeaturePos <- which(colnames(seurat_integrated@meta.data) == splitFeature)

## Workaround: subset does not allow for programmatical subsetting, therefore the column name cannot be provided as variable
## To solve this, I create a fake feature which is a copy of the selected one, but with a fixed name instead
seurat_integrated$myFeature <- seurat_integrated@meta.data[, splitFeaturePos]
seurat_integrated$myContrast <- seurat_integrated@meta.data[, contrastFeaturePos]
Idents(seurat_integrated) <- seurat_integrated$myContrast
message("Subsetting the seurat object by feature " %&% splitFeature)
allSplitFeatures <- unique(seurat_integrated$myFeature)
allContrastFeatures <- unique(seurat_integrated$myContrast)

allOutNames <- NULL
message("Beginning DE. Depending on data size each contrast can take several minutes...")
for (ii in seq_len(length(allSplitFeatures))) {
  mymarkers_pairwise <- NULL
  mymarkers_all <- NULL
  message("Working on " %&% allSplitFeatures[ii] %&% "\tFeature " %&% ii %&% " of " %&% length(allSplitFeatures))
  seurat_subset <- subset(seurat_integrated, subset = myFeature == allSplitFeatures[ii])
  kk <- 0

  for (jj in 1:(length(allContrastFeatures) - 1)) {
    for (kk in (jj + 1):length(allContrastFeatures)) {
      message("Contrast " %&% allContrastFeatures[jj] %&% " vs " %&% allContrastFeatures[kk])
      mymarkers_pairwise <- FindMarkers(object = seurat_subset, ident.1 = allContrastFeatures[jj], ident.2 = allContrastFeatures[kk], assay = myassay, logfc.threshold = 0, min.cells.group = 1, min.cells.feature = 1, min.pct = 0, only.pos = FALSE, test.use = test.use, latent.vars = latent.vars)
      outname <- allSplitFeatures[ii] %&% "_" %&% myassay %&% "-" %&% allContrastFeatures[jj] %&% "_vs_" %&% allContrastFeatures[kk]
      allOutNames <- c(allOutNames, outname)
      saveRDS(mymarkers_pairwise, outname)
      message("Plotting the boxplot for contrast " %&% allContrastFeatures[jj] %&% " vs " %&% allContrastFeatures[kk])
      mycounts <- as.data.frame(GetAssayData(object=seurat_subset, slot="counts", assay=myassay))
      mytopgenes <- head(mymarkers_pairwise[order(mymarkers_pairwise$p_val_adj),], 40)
      mycounts <- mycounts[which(rownames(mycounts) %in% rownames(mytopgenes)),]
      mycounts$gene = rownames(mycounts)
      mycounts_melt <- melt(mycounts)
      tmp_condi <- data.frame(condi = seurat_subset$myContrast, variable = names(seurat_subset$myContrast))
      mycounts_melt <- left_join(mycounts_melt, tmp_condi, by="variable")
      mycounts_melt <- mycounts_melt[which(mycounts_melt$condi %in% allContrastFeatures[c(jj,kk)]),]
      p1 <- ggplot(mycounts_melt, aes(x = condi, y = value, fill = condi)) +
        geom_boxplot(width = 0.4) +
        xlab(contrastFeature) +
        ylab("") +
        facet_wrap(~gene, nrow = 5, scales = "free_y") +
        theme(axis.text.x = element_blank()) +
        scale_fill_discrete(name = contrastFeature) +
        ggtitle(allSplitFeatures[ii] %&% ": " %&% contrastFeature)
      ggsave(outdir %&% allSplitFeatures[ii] %&% "__" %&% contrastFeature %&% "__" %&% paste(allContrastFeatures[c(jj,kk)], collapse="-vs-")  %&% "__de_boxplot.png",
        p1,
        width = 30, height = 15, units = "cm"
      )
      ###
    }
  }
  for (jj in seq_len(length(allContrastFeatures))) {
    message("Contrast " %&% allContrastFeatures[jj] %&% " vs all")
    mymarkers_all <- FindMarkers(object = seurat_subset, ident.1 = allContrastFeatures[jj], assay = myassay, logfc.threshold = 0, min.cells.group = 1, min.cells.feature = 1, min.pct = 0, only.pos = FALSE, test.use = test.use, latent.vars = latent.vars)
    allOutNames <- c(allOutNames, outname)
    outname <- allSplitFeatures[ii] %&% "_" %&% myassay %&% "-" %&% allContrastFeatures[jj] %&% "_vs_all"
    saveRDS(mymarkers_all, outname)
  }
}

message("Filtering results...")
for (i in seq_len(length(allOutNames))) {
  mydetable <- readRDS(allOutNames[i])
  ## If the loaded table is empty, nothing was found, which should not happen. Throw a warning
  if (nrow(mydetable) == 0) {
    warning("WARN: The DE table " %&% allOutNames[i] %&% " is empty. Please check the contrast and your data for inconsistencies as this table should include all unfiltered results")
  }
  mydetable <- mydetable[which(mydetable$p_val_adj < myfdr), ]
  mydetable <- mydetable[which(abs(mydetable$avg_log2FC) > myfc), ]
  if (nrow(mydetable) == 0) {
    message("The selected threshold lead to 0 selected markers. Please inspect the complete table in the RDS object " %&% allOutNames[i] %&% " to understand why")
  } else {
    message("Found " %&% nrow(mydetable) %&% " significant and strong markers for contrast " %&% allOutNames[i])
  }
  write.table(mydetable, paste0(allOutNames[i], "_filtered.tsv"), sep = "\t", quote = F)
}

message("Plotting images:")
for (i in seq_len(length(allOutNames))) {
  mydetable <- readRDS(allOutNames[i])
  mydetable$status <- "not"
  mydetable$status[(abs(mydetable$avg_log2FC) > myfc) & mydetable$p_val_adj < myfdr] <- "signif&strong"
  mydetable$status[(abs(mydetable$avg_log2FC) < myfc) & mydetable$p_val_adj < myfdr] <- "signif"
  mydetable$status[(abs(mydetable$avg_log2FC) > myfc) & mydetable$p_val_adj > myfdr] <- "strong"

  message("Plotting the volcano plot for " %&% allOutNames[i])
  p <- ggplot(data = mydetable, aes(x = avg_log2FC, y = -log10(p_val_adj), col = status)) +
    ggtitle(allOutNames[i]) +
    geom_point() +
    theme_minimal() +
    geom_hline(yintercept = -log10(myfdr)) +
    geom_vline(xintercept = c(myfc, -myfc)) +
    scale_color_manual(values = c("black", "blue", "red", "darkgreen"))
  ggsave(paste0(allOutNames[i], "_volcano.png"), plot = p)

  message("Plotting a zoomed-in version of the volcano plot for " %&% allOutNames[i])
  p <- ggplot(data = mydetable, aes(x = avg_log2FC, y = -log10(p_val_adj), col = status)) +
    ggtitle(allOutNames[i]) +
    geom_point() +
    theme_minimal() +
    geom_hline(yintercept = -log10(myfdr)) +
    geom_vline(xintercept = c(myfc, -myfc)) +
    scale_color_manual(values = c("black", "blue", "red", "darkgreen")) +
    coord_cartesian(xlim = c(-10, 10))
  ggsave(paste0(allOutNames[i], "_volcano_zoom.png"), plot = p)
}

  message("Plotting the Pvalue density for " %&% allOutNames[i])
  p <- ggplot(data = mydetable, aes(p_val, ..density.., fill="grey"))+
          geom_histogram(bins = 100) +
          geom_hline(yintercept = 1, col = "black") +
          xlab("P-value") +
          ylab("Density") +
          scale_fill_brewer(palette = "Dark2") +
          theme(legend.position = "none") +
          scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
          coord_cartesian(ylim = c(0, 5)) +
          ggtitle(allOutNames[i])
  ggsave(paste0(allOutNames[i], "_pvd.png")

  #TODO: complete the boxplots
  message("Plotting the boxplots for the top 40 DE genes")
  p <- ggplot(t_plot, aes(x = condi, y = value, fill = condi)) +
              geom_boxplot(width = 0.4) +
              xlab(this_comp) +
              ylab("") +
              facet_wrap(~gene, nrow = 5, scales = "free_y") +
              theme(axis.text.x = element_blank()) +
              scale_fill_manual(values = this_annot_colors, name = this_comp) +
              ggtitle(this_ct %&% ": " %&% this_comp)
            if (length(levels(t_plot$gene)) < 15) {
              p1 <- p1 + geom_jitter(width = 0.2, color = "grey50", size = 2)
            } else {
              p1 <- p1 + geom_jitter(width = 0.2, color = "grey50", size = 0.5)
            }
  ggsave(paste0(allOutNames[i], "_boxplots.png")
