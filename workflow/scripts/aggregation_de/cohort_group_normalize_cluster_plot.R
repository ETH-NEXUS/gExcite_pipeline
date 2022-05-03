#!/usr/bin/env Rscript
#################################################
## File name: cohort_group_normalize_cluster_plot.R
## Author: Matteo Carrara
## Date created: Oct 2021
## R Version: 4.1.0
##################################################

## GENERAL:
## This script is the second part of a 3-part processing of single-cell data. After the aggregation of the cohort samples performed with
## the previous script, the data can be normalized cluster-by-cluster.
## Several plots are produced alongside the actual output to provide insight on the data produced

### ADDITIONAL IMPORTANT INFORMATION:
## This script takes as input the results from the aggregation script and assumes the input follows the format and naming conventions
## introduced by the aggregation script. Please see the documentation below for detailed information on input and output
## All options are required and must have a value set. In the future there might be a migration away from optparse and to a package allowing to
## automatically check if a required option has a value or not

## INPUT:
## - The folder where the results of the aggregation step are stored
## - The grouping feature used to aggregate the samples
## - The directory in which the output must be stored
## - A custom name that will be used for some output files. It is suggested to use the cohort name
## - The minimum number of cells in a cluster below which the cluster is filtered out
## - The maximum number of highly-variable genes to use for the analysis
## - The variables to be used for testing. Up to a maximum of 7
## - The file containing the genes involved in the pathways of interest. GMT format or RDS object (Optional)
## - The file containing the marker genes for each celltype. GMX format (Optional)
## - The file containing the celltypes. Tab-delimited table with 'Major' and 'Subtype' celltypes as columns (Optional)
## - The file containing the clinical data. Tab-delimited table. The sample identifiers must be identical to those available in the aggregated RDS files (colnames)
## - The name of the column of the clinical data table containing the sample identifiers
## - The file containing the aggregate cell counts for each feature. This file is produced by the aggregation script

## OUTPUT:
## - For each element of the grouping feature (e.g. for each celltype if grouping is done by celltype)
##   - RDS object with the log-transformed normalized counts
##   - RDS object with the GSVA enrichment results
##   - Heatmap of celltype-specific genes (Optional)
##   - Heatmap of highly variable genes
##   - Expression scatterplot of the highly variable genes
##   - UMAP of the highly variable genes
##   - Violin plot of rlog-normalized data
##   - Table of GSVA enrichment distances
##   - Heatmap of GSVA results for the full cohort (Optional)

## ADDITIONAL INFORMATION:
##  - The sample identifiers are the only element that keeps track of the data throughout the analysis. It is very important to keep the identifiers
##    consistent between input files. The two major elements that make use of the sample identifiers are the aggregation RDS files and the
##    clinical data table
##  - celltype identifiers must be consistent in all input data. This means that the celltypes available in the aggregation filenames must be consistent with the celltypes listed in in "celltype_maker_genes" and "celltype_config"

## PATHWAYS file format
##    The pathways file should be in a standard gmt format, one pathway per line, followed by the genes in the ppathways, separated by a tab.. An RDS can be provided instead, but the object must still contain the gmt-formatted table as it is loaded in R

## CELLTYPE_MARKER_GENES file format
##    The file should be in a standard gmx format, on celltype each column, celltype name as column name and the rest of the column listing all genes part of the celltype

## CELLTYPE_CONFIG file format
##    The file should be a tab-delimited table with two columns, one listing the major celltypes and one listing all subtypes of that major celltype, separated by a comma. The column names must be "Major" and "Subtype". The names of the major celltypes listed here should be identical to the names of the major celltypes throughout the analysis

suppressPackageStartupMessages({
  library(DESeq2)
  library(reshape2)
  library(ggplot2)
  library(ggrepel)
  library(scran)
  library(dbscan)
  library(dynamicTreeCut)
  library(RColorBrewer)
  library(pheatmap)
  library(GSVA)
  library(zeallot)
  library(optparse)
})

option_list <- list(
  make_option("--aggregation_folder", type = "character", help = "Folder containing the results of the aggregation function"),
  make_option("--groupvar", type = "character", help = "Feature used to aggregate the samples. It will be used to retrieve the aggregation files."),
  make_option("--outdir", type = "character", help = "Full path to output directory."),
  make_option("--outName", type = "character", help = "Prefix name of output files."),
  make_option("--min_sample_size", type = "integer", default = 5, help = "The minimum sample size for the sample to be considered in the analysis"),
  make_option("--max_hvg_length", type = "integer", default = 1000, help = "Maximum number of high variable genes (hvg) to consider during the analysis"),
  make_option("--test_vars", type = "character", help = "Array containing all variables to test. The maximum amount of concurrent variables is 7. Additional variables will be dropped"),
  make_option("--pathways", type = "character", default = NULL, help = "File containing the gmt-formatted set of pathways with gene lists. Alternatively, an RDS containing the same information can be provided"),
  # make.option("--condition_marker_genes", type = "character", default = NULL,  help = "File containing a list of markers for the condition. One marker per line"),
  make_option("--celltype_marker_genes", type = "character", default = NULL, help = "File containing the gmx-formatted set of markers for each celltype"),
  make_option("--celltype_config", type = "character", default = NULL, help = "File containing the celltype configuration table. 2 columns called 'Major' and 'Subtype' respectively"),
  make_option("--clinical_data", type = "character", help = "File containing the clinical data. One column must contain the same identifiers used in the RDS files"),
  make_option("--clinical_data_id", type = "character", help = "Name of the column containing the identifiers in the clinical data table"),
  make_option("--meta_counts", type = "character", help = "File containing the aggregated counts for each groupvar")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

"%&%" <- function(a, b) paste(a, b, sep = "")

# WARNING: The clinid must be exactly the same as what comes as columns in the RDS from the aggregation!!!!!!!!!

#####################
## Initializing vars
#####################

message("Setting up the input")
aggregation_folder <- opt$aggregation_folder
groupvar <- opt$groupvar
outdir <- opt$outdir
outName <- opt$outName
min_sample_size <- opt$min_sample_size
max_hvg_length <- opt$max_hvg_length
test_vars <- opt$test_vars
gsfn <- opt$pathways
ctfn <- opt$celltype_marker_genes
ct_configfn <- opt$celltype_config
clinfn <- opt$clinical_data
clinid <- opt$clinical_data_id
metafn <- opt$meta_counts

rdsfn <- grep("RDS", grep(groupvar, list.files(aggregation_folder), value = TRUE), value = TRUE)
rdsfn <- aggregation_folder %&% rdsfn

gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

palette_set <- data.frame(
  pal_name = c(
    "Paired", "Dark2", "Set1", "Accent",
    "Set3", "Dark2", "Paired", "Accent"
  ),
  pal_length = c(12, 8, 8, 8, 12, 8, 12, 8)
)

########################
## Function declaration
########################

# Create the annotation colors and names for the plots
f_make_annot <- function(mytest_vars, myd_clin, myclinid) {
  annot <- as.data.frame(myd_clin[, mytest_vars])
  rownames(annot) <- myd_clin[, myclinid]
  for (mm in seq_len(ncol(annot))) {
    annot[, mm] <- factor(annot[, mm])
  }
  annot_colors <- vector("list", length(mytest_vars))
  names(annot_colors) <- mytest_vars
  for (nn in seq_along(annot_colors)) {
    this_pal <- brewer.pal(
      palette_set$pal_length[nn],
      palette_set$pal_name[nn]
    )[seq_along(levels(annot[, nn]))]
    annot_colors[[mytest_vars[nn]]] <- this_pal
    names(annot_colors[[nn]]) <- levels(annot[, nn])
  }
  # HERE IT makes a difference between what was the tupro id and the patient id. This is not always the case and seems specific
  annot$patient <- factor(myd_clin[, myclinid])
  annot_colors$patient <- gg_color_hue(length(levels(annot$patient)))
  names(annot_colors$patient) <- levels(annot$patient)
  annot <- annot[, c("patient", colnames(annot)[seq_len(ncol(annot) - 1)])]
  return(list(annot, annot_colors))
}

#################
## Main code
#################

# load clinical annotation
message("Loading the clinical annotation")
d_clin <- read.table(clinfn, head = TRUE, sep = "\t", check.names = FALSE)
names(d_clin) <- gsub(" |-", "_", names(d_clin))
d_clin[d_clin == ""] <- NA

test_vars <- read.delim(test_vars, header = FALSE)
test_vars <- test_vars[, ncol(test_vars)]
test_vars <- test_vars[seq_len(min(7, length(test_vars)))]

c(annot, annot_colors) %<-% f_make_annot(test_vars, d_clin, clinid)

# rlog-transform per group, plots, clustering
message("Loading the cell counts file")
meta_df <- read.table(metafn, head = TRUE, sep = "\t")

meta_df$sample <- gsub("(^.*)-.*_(.*)", "\\1_\\2", meta_df$sample)
meta_df$sample <- gsub("(^.*)-.*", "\\1", meta_df$sample)
str(meta_df)

for (ii in seq_along(rdsfn)) {
  message(paste("Running the normalization on file", rdsfn[[ii]], "\t", ii, "of", length(rdsfn)))
  message("Loading the input aggregation file")
  dat <- readRDS(rdsfn[ii])
  dat <- as.matrix(dat, drop = FALSE)
  col_names <- gsub("(^.*)-.*_(.*)", "\\1_\\2", colnames(dat))
  col_names <- gsub("(^.*)-.*", "\\1", col_names)
  id_keep <- !duplicated(col_names)
  dat <- dat[, id_keep, drop = FALSE]
  colnames(dat) <- col_names[id_keep]
  this_ct <- unlist(strsplit(gsub(aggregation_folder, "", rdsfn[ii]), split = "_"))[1]
  this_meta <- meta_df[which(meta_df$cluster == this_ct), ]
  id_keep <- !duplicated(this_meta$sample)
  this_meta <- this_meta[id_keep, , drop = FALSE]
  rownames(this_meta) <- this_meta$sample
  dat <- dat[, this_meta$sample, drop = FALSE]
  # omit samples w/o clinical info
  sampleID <- gsub("(.*)_.*", "\\1", colnames(dat))
  idx <- match(sampleID, rownames(annot))
  idx <- idx[!is.na(idx)]
  this_samples <- sampleID %in% rownames(annot)[idx]
  dat <- dat[, which(this_samples), drop = FALSE]
  this_meta <- this_meta[colnames(dat), , drop = FALSE]
  if (ncol(dat) > min_sample_size) {
    # normalize with rlog / vst
    message("Normalizing the data")
    logdat <- varianceStabilizingTransformation(dat)
    d_plot <- melt(logdat, varnames = c("gene", "sample"))
    d_plot$type <- "rlog"
    tmp <- melt(dat, varnames = c("gene", "sample"))
    tmp$type <- "log_raw"
    tmp$value <- log2(tmp$value + 1)
    d_plot <- rbind(d_plot, tmp)
    d_plot$type <- factor(d_plot$type, levels = c("log_raw", "rlog"))
    message("Plotting the violin plot")
    rlog_violin_plot <- ggplot(d_plot, aes(x = sample, y = value, fill = sample)) +
      ggtitle(this_ct) +
      geom_violin() +
      facet_wrap(~type, nrow = 1) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none"
      )
    ggsave(outdir %&% this_ct %&% "__" %&% groupvar %&% "__rlog_violin.png",
      width = 18, height = 8, units = "cm", plot = rlog_violin_plot
    )

    # hvg detection
    # var.vs.mean is not always constant!
    message("High-variable genes detection and expression plotting")
    meanvartrend <- modelGeneVar(logdat)
    hvg <- getTopHVGs(meanvartrend, fdr.threshold = 0.05)
    gene_desc <- meanvartrend
    gene_desc$symbol <- rownames(meanvartrend)
    gene_desc$hvg <- (gene_desc$symbol %in% hvg)
    hvg_plot <- ggplot(as.data.frame(gene_desc), aes(x = mean, y = total)) +
      geom_point(aes(color = hvg)) +
      geom_smooth() +
      geom_line(aes(x = mean, y = tech), color = "darkgreen", size = 2)
    ggsave(outdir %&% this_ct %&% "__" %&% groupvar %&% "__hvg_mean-total.png",
      width = 18, height = 8, units = "cm", plot = hvg_plot
    )
    # pick max number of hvg's
    dtest <- gene_desc[hvg, ]
    dtest <- dtest[order(dtest$bio, decreasing = TRUE), ]
    dtest <- dtest[seq_len(min(nrow(dtest), max_hvg_length)), ]
    gene_desc[setdiff(rownames(gene_desc), rownames(dtest)), ]$hvg <- FALSE
    hvg <- rownames(gene_desc)[which(gene_desc$hvg)]
    # umap based on hvg
    message("Calculating and plotting the UMAP based on the hvgs")
    rumap <- umap(t(logdat[hvg, ]), n_neighbors = 5, scale = TRUE, ret_nn = TRUE)
    embed <- as.data.frame(rumap$embedding)
    rownames(embed) <- colnames(logdat)
    names(embed) <- c("umap_1", "umap_2")
    dbc <- hdbscan(embed, minPts = 3)
    embed$dbscan <- factor(dbc$cluster)
    this_meta <- cbind(this_meta, embed[this_meta$sample, ])
    ggp <- ggplot(embed, aes(x = umap_1, y = umap_2, color = dbscan, label = rownames(embed))) +
      theme_bw() +
      geom_label_repel(show.legend = FALSE, max.overlaps = 20) +
      geom_point() +
      ggtitle(this_ct) +
      guides(color = guide_legend(override.aes = list(size = 3))) 
    if (length(unique(embed$dbscan)) > 8) ggp <- ggp + theme(legend.position = "none")
    ggsave(outdir %&% this_ct %&% "__" %&% groupvar %&% "__hvg_umap.png",
      width = 16, height = 16, units = "cm", plot = ggp
    )

    # heatmap of hvgs
    message("Calculating and plotting the heatmap of the hvgs")
    t_labels <- rownames(gene_desc)[order(gene_desc$bio, decreasing = TRUE)][seq_len(min(length(hvg), 20))]
    d_plot <- logdat[hvg, ]
    rownames(d_plot)[!rownames(d_plot) %in% t_labels] <- ""
    # annotation
    sampleID <- gsub("(.*)_.*", "\\1", colnames(d_plot))
    idx <- match(sampleID, rownames(annot))
    idx <- idx[!is.na(idx)]
    this_samples <- sampleID %in% rownames(annot)[idx]
    d_plot <- d_plot[, which(this_samples)]
    this_annot <- annot[idx, ]
    rownames(this_annot) <- colnames(d_plot)
    col_labels <- d_clin[, clinid][match(rownames(annot)[idx], d_clin[, clinid])]
    hmcol <- rev(colorRampPalette(brewer.pal(8, "RdBu"))(127))


    pat_lab <- names(which(table(this_annot$patient) > 1))
    if (length(pat_lab) > 0 & groupvar != "ct_cl") {
      this_annot$patient[!this_annot$patient %in% pat_lab] <- NA
    }
    if (groupvar == "ct_cl") {
      col_labels <- rep("", length(col_labels))
    }
    this_annot$patient <- factor(this_annot$patient)
    annot_colors$patient <- gg_color_hue(length(levels(this_annot$patient)))
    names(annot_colors$patient) <- levels(this_annot$patient)

    hc <- hclust(dist(t(d_plot)), method = "ward.D2")
    ccl <- cutreeHybrid(hc, as.matrix(dist(t(d_plot))), minClusterSize = 3)
    this_meta$hvg_cl <- NA
    this_meta[colnames(d_plot), "hvg_cl"] <- ccl$labels
    this_meta$hvg_cl <- factor(this_meta$hvg_cl)
    phm <- pheatmap(d_plot,
      scale = "none", clustering_method = "ward.D2",
      cluster_cols = hc, cutree_cols = length(levels(this_meta$hvg_cl)),
      annotation_col = this_annot, annotation_colors = annot_colors,
      show_rownames = TRUE, show_colnames = TRUE, fontsize = 8,
      labels_col = col_labels,
      color = hmcol, main = this_ct %&% " HVGs"
    )
    ggsave(outdir %&% this_ct %&% "__" %&% groupvar %&% "__hvg_hm.png",
      phm$gtable,
      width = 30, height = 30, units = "cm"
    )

    # heatmap of celltype specific genes, ONLY if the file is provided
    if (!is.null(ctfn)) {
      # load celltype marker gene list
      message("Loading the celltype marker gene list")
      ct_genes <- read.table(ctfn, head = TRUE, fill = TRUE, sep = "\t")
      ct_genes <- apply(ct_genes, 2, function(x) x[x != ""])
      ct_config <- read.table(ct_configfn, head = TRUE, sep = "\t")
      ct_genes <- ct_genes[ct_config$Major]
      ## This is very specific: it takes composed major celltypes names, e.g. "T.cells_melanoma_Tirosh16" and reduces them to just the celltype, e.g. "T.cells". This is not necessary if consistency between file content is kept
      names(ct_genes) <- sapply(names(ct_genes), function(x) strsplit(x, "_")[[1]][1])

      message("Calculating and plotting the heatmap of celltype-specific genes")

      idx <- match(ct_genes[[this_ct]], rownames(logdat))
      d_plot <- logdat[idx[!is.na(idx)], which(this_samples)]
      # omit genes not changing across samples
      highvar <- which(rowVars(d_plot) > 1e-6)
      hc <- hclust(dist(t(d_plot)), method = "ward.D2")
      ccl <- cutreeHybrid(hc, as.matrix(dist(t(d_plot))), minClusterSize = 3)
      this_meta$ctg_cl <- NA
      this_meta[colnames(d_plot), "ctg_cl"] <- ccl$labels
      this_meta$ctg_cl <- factor(this_meta$ctg_cl)
      phm <- pheatmap(d_plot[highvar, ],
        scale = "none", clustering_method = "ward.D2",
        cluster_cols = hc, cutree_cols = length(levels(this_meta$ctg_cl)),
        annotation_col = this_annot, annotation_colors = annot_colors,
        fontsize = 8, show_rownames = TRUE, show_colnames = TRUE,
        labels_col = col_labels,
        color = hmcol, main = this_ct %&% " Marker Genes"
      )
      ggsave(outdir %&% this_ct %&% "__" %&% groupvar %&% "__ctg_hm.png",
        phm$gtable,
        width = 50, height = 35, units = "cm"
      )
    } else {
      message("Celltype marker genes not provided. Skipping related heatmap")
    }

    # gsva, ONLY if the file is provided
    if (!is.null(gsfn)) {
      # Load the gene lists. They can either be RDS or gmt. This code will handle both
      message("Loading the pathway file")
      is_rds <- endsWith(gsfn, ".RDS")
      if (is_rds) {
        gsidl <- readRDS(gsfn)
      }
      is_gmt <- endsWith(gsfn, ".gmt")
      if (is_gmt) {
        # load geneset list
        gsidl <- readLines(gsfn)
        gsidl <- lapply(gsidl, function(x) strsplit(x, "\\\t")[[1]])
        names(gsidl) <- sapply(gsidl, function(x) x[1])
        # remove gene set names from gene list
        gsidl <- sapply(gsidl, function(x) x[-1])
        # remove gene set description from gene list
        gsidl <- sapply(gsidl, function(x) x[-1])
        names(gsidl) <- gsub("HALLMARK_", "", names(gsidl))
      }

      message("Calculating and plotting gsva")
      rgsva <- gsva(logdat, gsidl, method = "gsva")
      d_plot <- rgsva[, which(this_samples)]
      clone_dist <- dist(t(d_plot))
      colclust <- hclust(clone_dist, method = "ward.D2")
      ccl <- cutreeHybrid(colclust, as.matrix(dist(t(d_plot))), minClusterSize = 3)
      ccl$labels <- factor("ph" %&% (ccl$labels))
      rowclust <- hclust(dist(d_plot), method = "ward.D2")
      rcl <- cutreeHybrid(rowclust, as.matrix(dist(d_plot)), minClusterSize = 3)
      rcl$labels <- factor("si" %&% (rcl$labels))
      phm <- pheatmap(d_plot,
        scale = "none", color = hmcol, fontsize = 8, fontsize_row = 6,
        cluster_cols = colclust, cutree_cols = length(levels(ccl$labels)),
        cluster_rows = rowclust, cutree_rows = length(levels(rcl$labels)),
        annotation_colors = annot_colors, annotation_col = this_annot,
        show_rownames = T, show_colnames = T, labels_col = col_labels,
        main = this_ct %&% ": " %&% outName %&% " geneset"
      )
      ggsave(outdir %&% this_ct %&% "__" %&% groupvar %&% "__" %&% outName %&% "__" %&% mypwname %&% "__gsva_hm.png",
        phm$gtable,
        width = 40, height = 50, units = "cm"
      )
      clone_dist <- as.matrix(clone_dist)
      rownames(clone_dist) <- colnames(clone_dist) <- colnames(d_plot)
      mypwname <- unlist(lapply(strsplit(gsfn, "/"), function(x) x[length(x)]))
      mypwname <- gsub("\\.gmt", "", mypwname)
      write.table(clone_dist, outdir %&% this_ct %&% "__" %&% groupvar %&% "__" %&% outName %&%
        "__" %&% mypwname %&% "__gsva_distm.txt", col.names = NA, row.names = TRUE, sep = "\t")
    } else {
      message("Pathway genes not provided. Skipping related heatmap")
    }

    # build SummarizedExperiment object and save.
    message("Building and saving the RDS output")
    this_meta <- cbind(this_meta, d_clin[match(this_meta$sample, d_clin[, clinid]), ])
    tmp <- strsplit(clinfn, split = "/")[[1]]
    clin_fn <- tmp[length(tmp)]
    se <- SummarizedExperiment(
      assays = SimpleList(normcounts = logdat),
      rowData = gene_desc,
      colData = DataFrame(this_meta),
      metadata = list(
        groupvar = groupvar,
        groupval = this_ct,
        clin_fn = clin_fn,
        umap_nn = rumap$nn
      )
    )
    saveRDS(se, outdir %&% this_ct %&% "__" %&% groupvar %&% "__se.RDS")

    if (!is.null(gsfn)) {
      se <- SummarizedExperiment(
        assays = SimpleList(gsva = rgsva),
        colData = DataFrame(this_meta),
        metadata = list(geneset = gsidl)
      )
      saveRDS(se, outdir %&% this_ct %&% "__" %&% groupvar %&% "__" %&% outName %&% "__gsva_se.RDS")
    }
  }
}

message("Run complete!")
