#!/usr/bin/env Rscript
#################################################
## File name: ....
## Author: Matteo Carrara
## Date created: Aug 2021
## R Version: 4.1.0
##################################################

## TESTED ON THE FOLLOWING VERSIONS OF THE DEPENDENCIES:
## multcomp 1.4-18
## DESeq2 1.34.0
## reshape2 1.4.4
## ggplot2 3.3.5
## ggrepel 0.9.1
## scran 1.22.0
## dynamicTreeCut 1.63-1
## RColorBrewer 1.1-2
## pheatmap 1.0.12
## GSVA 1.42.0
## UpSetR 1.4.0
## zeallot 0.1.0
## optparse 1.7.1

## GENERAL:
## This script is the third part of a 3-part processing of single-cell data. After the aggregation and normalization of the cohort samples performed with
## the previous scripts, the data can be analyzed for differential expression cluster-by-cluster to identify DE genes in the different test variable chosen.
## Several plots are produced alongside the actual output to provide insight on the data produced

### ADDITIONAL IMPORTANT INFORMATION:
## This script takes as input the results from the ... script and assumes the input follows the format and naming conventions
## introduced by the previous scripts. Please see the documentation below for detailed information on input and output
## All options are required and must have a value set. In the future there might be a migration away from optparse and to a package allowing to
## automatically check if a required option has a value or not

## INPUT:
## - The folder where the results of the aggregation step are stored. Please provide the main folder
## - The grouping feature used to aggregate the samples
## - The directory in which the output must be stored
## - The minimum number of cells in a cluster below which the cluster is filtered out
## - The variables to be used for testing. Text table with 1 contrast per line. Up to a maximum of 7 lines. Please refer to section "TEST VARIABLES TABLE for details about the structure
## - The file containing the clinical data. Tab-delimited table. The sample identifiers must be identical to those available in the aggregated RDS files (colnames)
## - The name of the column of the clinical data table containing the sample identifiers
## - The FDR cutoff to use to define a gene as significant
## - The fold-change cutoff to use to define a gene as strong

## OUTPUT:
## - For each element of the grouping feature (e.g. for each celltype if grouping is done by celltype)
##   AND for each test variable from the clinical data:
##   - Boxplot of gene expression of the most significant genes (if any). "Most significant" means that they pass both the expression and the FDR cutoffs
##   - Heatmap of gene expression of the most significant genes (if any)
##   - Vennplot (UpSet R) showing in which contrasts we can find the different most significant genes (if any)
##   - PVD plot (Pvalue density barplot)
##   - Volcano plot
##   - Dispersion estimate plot

## REMARK on sample identifiers
## The sample identifiers are the only element that keeps track of the data throughout the analysis. It is very important to keep the identifiers
## consistent between input files. The two major elements that make use of the sample identifiers are the RDS files and the
## clinical data table

## TEST VARIABLE TABLE
## It must be tab delimited, with no row names and no column names. Every line contains one specific contrast, one element per cell. The last column must always contain the main variable of interest. All previous cells contain possible covariates.
## If the contrasts in the table have different numbers of covariates, empty covariate cells are acceptable and will be ignored.
## All variable names must match column names available in the clinical data file.
## The script prints on screen the final model. The contrast matrix is built based on the main variable. One column per value of each covariate is added to the contrast matrix with all contrasts set to 0.
## Please refer to the file "test_vars_example.txt" for an example of the table structure.
##

suppressPackageStartupMessages({
  library(multcomp)
  library(DESeq2)
  library(reshape2)
  library(ggplot2)
  library(ggrepel)
  library(scran)
  library(dynamicTreeCut)
  library(RColorBrewer)
  library(pheatmap)
  library(GSVA)
  library(UpSetR)
  library(zeallot)
  library(optparse)
})

option_list <- list(
  make_option("--aggregation_folder", type = "character", help = "Folder containing the results of the aggregation function"),
  make_option("--groupvar", type = "character", help = "Feature used to aggregate the samples. It will be used to retrieve the aggregation files."),
  make_option("--outdir", type = "character", help = "Full path to output directory."),
  make_option("--min_sample_size", type = "integer", default = 5, help = "The minimum sample size for the sample to be considered in the analysis"),
  make_option("--test_vars", type = "character", help = "Filename of a text file containing all variables to test. One variable name per line. The maximum amount of concurrent variables is 7. Additional variables will be dropped"),
  make_option("--clinical_data", type = "character", help = "File containing the clinical data. One column must contain the same identifiers used in the RDS files"),
  make_option("--clinical_data_id", type = "character", help = "Name of the column containing the identifiers in the clinical data table"),
  # make_option("--meta_counts", type = "character", help = "File containing the aggregated counts for each groupvar"),
  # make_option("--verbose", type = "integer", default = 1, help = "The level of verbosity of the script"),
  make_option("--fdr.cut", type = "double", default = 0.05, help = "The FDR cutoff for the DE analysis"),
  make_option("--fc.cut", type = "double", default = 1.5, help = "The fold-change cutoff for the DE analysis")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

"%&%" <- function(a, b) paste(a, b, sep = "")

# if (opt$verbose == 1) {
message("Setting up the input")
# }
aggregation_folder <- opt$aggregation_folder %&% "/"
groupvar <- opt$groupvar
outdir <- opt$outdir %&% "/"
min_sample_size <- opt$min_sample_size
test_vars_file <- opt$test_vars
test_vars <- read.delim(test_vars_file, header = FALSE, stringsAsFactors = FALSE)
clinfn <- opt$clinical_data
clinid <- opt$clinical_data_id
metafn <- aggregation_folder %&% "/aggregate." %&% groupvar %&% ".txt"
fdr.cut <- opt$fdr.cut
fc.cut <- opt$fc.cut

rdsfn <- grep((groupvar %&% ".RDS"), list.files(aggregation_folder), value = TRUE)
rdsfn <- aggregation_folder %&% rdsfn


theme_set(theme_bw(base_size = 10))

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


f_make_annot <- function(mytest_vars, myd_clin, myclinid) {
  annot <- as.data.frame(myd_clin[, mytest_vars])
  colnames(annot) <- mytest_vars
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
  annot$ID <- factor(myd_clin[, myclinid])
  annot_colors$ID <- gg_color_hue(length(levels(annot$ID)))
  names(annot_colors$ID) <- levels(annot$ID)
  annot <- annot[, c("ID", colnames(annot)[seq_len(ncol(annot) - 1)])]
  return(list(annot, annot_colors))
}

# load clinical annotation
message("Loading the clinical annotation")
d_clin <- read.table(clinfn, head = TRUE, sep = "\t", check.names = FALSE)
names(d_clin) <- gsub(" |-", "_", names(d_clin))
d_clin[d_clin == ""] <- NA
if (nrow(test_vars) > 7) {
  warning("There are more than 7 test variables. Working on the first 7")
}
test_vars <- as.data.frame(test_vars[seq_len(min(7, nrow(test_vars))), ])

# FIXME: This works only if there is one row in the test vars. It may be a good idea to require only one set of variables and repeat the entire script for multiple ones
c(annot, annot_colors) %<-% f_make_annot(unlist(test_vars), d_clin, clinid)


# sample inclusion matrix, 1 column per comparison, 1 line per sample:
message("Creating the sample inclusion matrix")
sample_incl <- matrix(FALSE, ncol = nrow(test_vars), nrow = nrow(d_clin))
rownames(sample_incl) <- d_clin[, clinid]
colnames(sample_incl) <- test_vars[, ncol(test_vars)]
for (ii in 1:(nrow(test_vars))) {
  tab <- table(d_clin[, test_vars[ii, ncol(test_vars)]])
  tab <- tab[tab > min_sample_size]
  idx <- which(names(tab) != "" & !is.na(names(tab)))
  samp_in <- d_clin[d_clin[, test_vars[ii, ncol(test_vars)]] %in% names(tab[idx]), clinid]
  sample_incl[samp_in, ii] <- TRUE
}


# DESeq2 per group, contrasts, plots,
message("Loading the cell counts file")
meta_df <- read.table(metafn, head = TRUE, sep = "\t")
meta_df$sample <- gsub("(^.*)-.*_(.*)", "\\1_\\2", meta_df$sample)
str(meta_df)
for (ii in seq_along(rdsfn)) {
  message(paste("Running DESeq2 on file", rdsfn[[ii]], "\t", ii, "of", length(rdsfn)))
  message("Loading the input RDS file")
  dat <- readRDS(rdsfn[ii])
  if (ncol(dat) >= 2 * min_sample_size) {
    dat <- as.matrix(dat, drop = FALSE)
    colnames(dat) <- gsub("(^.*)-.*_(.*)", "\\1_\\2", colnames(dat))
    this_ct <- unlist(strsplit(gsub(aggregation_folder, "", rdsfn[ii]), split = "_"))[1]
    this_meta <- meta_df[meta_df$cluster == this_ct, ]
    rownames(this_meta) <- this_meta$sample
    dat <- dat[, this_meta$sample, drop = FALSE]
    for (jj in seq_len(nrow(test_vars))) {
      message(paste("Working on test variable", test_vars[jj, ncol(test_vars)], "\t", jj, "of", nrow(test_vars)))
      # select samples for this comparison
      this_comp <- test_vars[jj, ncol(test_vars)]
      this_comp_full <- as.data.frame(test_vars[jj, ])
      if (length(which(this_comp_full == "")) != 0) {
        this_comp_full <- this_comp_full[-which(this_comp_full == "")]
      }
      this_samp <- match(names(which(sample_incl[, jj])), colnames(dat))
      this_samp <- this_samp[!is.na(this_samp)]
      this_dat <- dat[, this_samp]
      this_samp <- colnames(this_dat)
      idx <- match(this_samp, d_clin[, clinid])
      this_meta <- d_clin[idx, c(clinid, unlist(this_comp_full))]
      tab <- table(this_meta[, this_comp])
      tab <- tab[tab >= min_sample_size]
      for (aa in 2:ncol(this_meta)) {
        this_meta[, aa] <- factor(this_meta[, aa])
      }
      if (length(tab) >= 2) {
        message(paste("Test variable", test_vars[jj, ncol(test_vars)], "has 2 or more values above the minimum sample size threshold:\nPerforming DE analysis..."))
        # model formula for this comparison
        # dea_model <- formula("~" %&% this_comp)

        # Accept custom constrasts in input from the user. They must be columns from the clinical data table. At the moment, no sanity check is done about the structure of the contrast
        tmp <- paste(test_vars[jj, ])
        if (length(which(tmp == "")) != 0) {
          tmp <- tmp[-which(tmp == "")]
        }
        tmp <- paste(tmp, collapse = "+")
        dea_model <- formula("~" %&% tmp)
        message(paste("Based on the provided table, the model in position", jj, "is", paste(dea_model, collapse = "")))

        # DESeq object
        A <- rowMeans(this_dat)
        this_dat <- this_dat[A > 5, ]
        dds <- DESeqDataSetFromMatrix(
          countData = this_dat,
          colData = this_meta,
          design = dea_model
        )
        # Fit the model
        des <- DESeq(dds, betaPrior = T, quiet = TRUE)
        message("Plotting dispersion estimates")
        png(outdir %&% this_ct %&% "__" %&% this_comp %&% "__disp_ests.png", width = 30, height = 20, res = 600, units = "cm")
        plotDispEsts(des)
        dev.off()
        # contrast matrix for this comparison: all pairwise
        t_values <- NULL
        t_values <- levels(this_meta[, ncol(this_meta)])
        t_levels <- rep(1, length(t_values))
        names(t_levels) <- this_comp %&% t_values
        cmat <- contrMat(t_levels, type = "Tukey")
        cmat <- cbind("Intercept" = rep(0, nrow(cmat)), cmat)
        colnames(cmat) <- make.names(colnames(cmat))
        rownames(cmat) <- gsub(" - ", ".vs.", rownames(cmat))
        rownames(cmat) <- gsub(this_comp, "", rownames(cmat))
        rownames(cmat) <- make.names(rownames(cmat))
        rownames(cmat) <- gsub("\\.\\.", "\\.", rownames(cmat))
        colnames(cmat) <- gsub("\\.\\.", "\\.", colnames(cmat))
        rownames(cmat) <- gsub("\\.$", "", rownames(cmat))
        colnames(cmat) <- gsub("\\.$", "", colnames(cmat))

        if (ncol(this_comp_full) > 1) {
          for (bb in 1:(ncol(this_comp_full) - 1)) {
            tmp <- levels(this_meta[, which(colnames(this_meta) == this_comp_full[, bb])])
            for (cc in seq_len(length(tmp))) {
              cmat <- cbind(cmat, 0)
              colnames(cmat)[ncol(cmat)] <- paste0(this_comp_full[, bb], tmp[cc])
            }
          }
        }
        t_levels <- levels(this_meta[, this_comp])
        # extract contrasts
        dres <- data.frame(NULL)
        for (ll in seq_len(nrow(cmat))) {
          tres <- results(des, contrast = cmat[ll, ])
          tres <- as.data.frame(tres)
          tres$comparison <- rownames(cmat)[ll]
          tres$gene <- rownames(tres)
          rownames(tres) <- NULL
          dres <- rbind(dres, tres)
        }
        rm(tres)
        dres$padj <- p.adjust(dres$pvalue, method = "BH")
        dres$padj[dres$padj == 0] <- min(dres$padj[dres$padj > 0])
        dres$comparison <- factor(dres$comparison, levels = rownames(cmat))
        names(dres)[names(dres) == "log2FoldChange"] <- "logFC"
        # PVDs
        message("Plotting the PVDs")
        myplot <- ggplot(dres, aes(pvalue, ..density.., fill = "grey")) +
          geom_histogram(bins = 100) +
          geom_hline(yintercept = 1, col = "black") +
          xlab("P-value") +
          ylab("Density") +
          scale_fill_brewer(palette = "Dark2") +
          theme(legend.position = "none") +
          scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
          coord_cartesian(ylim = c(0, 5)) +
          facet_wrap(~comparison, ncol = 4) +
          ggtitle(this_ct %&% ": " %&% this_comp)
        ggsave(outdir %&% this_ct %&% "__" %&% this_comp %&% "__de_pvd.png",
          width = 30, height = 20, dpi = 600, units = "cm", plot = myplot
        )
        # Volcano plots
        message("Plotting the Volcano plot")
        dres$hit <- c("not", "signif")[as.numeric(dres$padj < fdr.cut) + 1]
        dres$hit[abs(dres$logFC) > log2(fc.cut) & dres$hit == "signif"] <- "signif & strong"
        dres$hit[abs(dres$logFC) > log2(fc.cut) & dres$hit == "not"] <- "strong"
        dres$hit <- factor(dres$hit, levels = c("not", "strong", "signif", "signif & strong"))
        cbPalette <- c(
          "#000000", "#009E73", "#56B4E9", "#D55E00",
          "#F0E442", "#0072B2", "#E69F00", "#CC79A7"
        )
        myvolc <- ggplot(dres[!is.na(dres$hit), ], aes(x = logFC, y = -log10(padj), col = hit)) +
          geom_point(shape = 16, size = 1) +
          geom_hline(yintercept = -log10(fdr.cut)) +
          geom_vline(xintercept = c(-1, 1) * log2(fc.cut)) +
          xlab(expression(paste("Effect size [", log[2](fold - change), "]"))) +
          ylab(expression(paste("Significance [", -log[10](padj), "]"))) +
          scale_color_manual(values = cbPalette, name = "DEA class", drop = FALSE) +
          theme(legend.justification = c(0, 1), legend.position = "top") +
          guides(color = guide_legend(override.aes = list(size = 3))) +
          ggtitle(this_ct %&% ": " %&% this_comp) +
          # geom_text(aes(label = ifelse(hit == "signif & strong", as.character(gene), "")), hjust = 0, vjust = 0) +
          facet_wrap(~comparison, ncol = 4) # , scales="free_y")
        ggsave(outdir %&% this_ct %&% "__" %&% this_comp %&% "__de_volcano.png",
          width = 30, height = 20, dpi = 600, units = "cm", plot = myvolc
        )
        # Top x most significant genes

        # Heatmap of top y most significant genes
        tops <- dres[!is.na(dres$gene) & !is.na(dres$pvalue), ]
        tops <- tops[order(tops$padj, decreasing = FALSE), ]
        tops <- tops[tops$hit == "signif & strong", ]
        saveRDS(tops, outdir %&% this_ct %&% "_" %&% this_comp %&% "_tops.RDS")
        sel <- unique(tops$gene)
        if (length(sel) > 0) {
          if (length(sel) <= 60) {
            tops <- sel
          } else {
            tmp <- sort(table(tops$comparison))
            choose_min <- tmp[1]
            sel <- as.character(NULL)
            while (length(sel) < 40) {
              for (i1 in seq_len(nrow(cmat))) {
                id_take <- which(tops$comparison == rownames(cmat)[i1])
                chosen <- setdiff(tops$gene[id_take], sel)
                id_rem <- seq_len(min(length(chosen), choose_min))
                sel <- c(sel, chosen[id_rem])
              }
              tops <- tops[which(!tops$gene %in% sel), ]
              tops$comparison <- factor(tops$comparison)
              tmp <- sort(table(tops$comparison))
              choose_min <- tmp[1]
            }
            tops <- unique(sel)[seq_len(40)]
          }
          if (length(tops) > 0) {
            message(paste("Test variable", test_vars[jj, ncol(test_vars)], "shows", length(tops), "strong and significant genes. We can continue"))
            # get row ids of all genes in tops in matrix mrld
            rld <- rlog(dds, blind = TRUE)
            mrld <- assay(rld)
            saveRDS(mrld, outdir %&% this_ct %&% "_" %&% this_comp %&% "_mrld.RDS")
            sel_x <- match(tops[seq_len(min(length(tops), 40))], rownames(mrld))
            sel_y <- match(colnames(this_dat), colnames(mrld))
            m_plot <- mrld[sel_x, sel_y, drop = F]
            # annotation
            idx <- match(colnames(mrld), d_clin[, clinid])
            idx <- idx[!is.na(idx)]

            # annotation
            sampleID <- colnames(m_plot)
            idx <- match(sampleID, rownames(annot))
            idx <- idx[!is.na(idx)]
            this_samples <- sampleID %in% rownames(annot)[idx]
            m_plot <- m_plot[, which(this_samples), drop = FALSE]
            this_annot <- annot[idx, ]
            rownames(this_annot) <- colnames(m_plot)
            col_labels <- d_clin[match(rownames(annot)[idx], d_clin[, clinid]), clinid]
            hmcol <- rev(colorRampPalette(brewer.pal(8, "RdBu"))(127))

            # boxplot of top significant genes per group
            message("Calculating and plotting the boxplot for the significant genes")
            t_plot <- melt(m_plot)
            names(t_plot)[1:2] <- c("gene", "sample")
            t_plot <- merge(t_plot, cbind(annot, sample = rownames(annot)))
            t_plot$condi <- t_plot[, this_comp]
            t_plot$gene <- factor(t_plot$gene, levels = sort(levels(t_plot$gene)))
            # Without removing references to other levels of the annotation from t_plot and annot_colors, the barplot
            # has additional useless elements in the legend
            t_plot[, this_comp] <- droplevels(t_plot[, this_comp])
            this_annot_colors <- annot_colors[[this_comp]]
            this_annot_colors <- this_annot_colors[which(names(this_annot_colors) %in% levels(t_plot[, this_comp]))]
            p1 <- ggplot(t_plot, aes(x = condi, y = value, fill = condi)) +
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
            ggsave(outdir %&% this_ct %&% "__" %&% this_comp %&% "__de_boxplot.png",
              p1,
              width = 30, height = 15, units = "cm"
            )

            # heatmap
            if (length(tops) > 1) {
              message("Test variable has more than 1 strong and significant gene. We can continue")
              message("Plotting the heatmap")
              phm <- pheatmap(t(m_plot),
                color = hmcol, scale = "column", fontsize = 12,
                cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "ward.D",
                show_colnames = TRUE, show_rownames = FALSE, annotation_names_col = TRUE,
                annotation_row = this_annot, annotation_colors = annot_colors,
                main = this_ct %&% ": " %&% this_comp
              )
              ggsave(outdir %&% this_ct %&% "__" %&% this_comp %&% "__de_hm.png", phm$gtable,
                width = 40, height = 30, units = "cm"
              )
              saveRDS(m_plot, outdir %&% this_ct %&% "_" %&% this_comp %&% "m_plot.RDS")
              # upset plot
              message("Plotting the upset plot")
              test_results <- dcast(dres[
                which(dres$hit == "signif & strong"),
                c("comparison", "gene")
              ],
              gene ~ comparison,
              fun = length
              )
              if (ncol(test_results) > 2) {
                p_upset <- upset(test_results[, -1],
                  order.by = "freq", nintersects = NA,
                  mainbar.y.label = this_ct %&% ": " %&% this_comp
                )
                png(outdir %&% this_ct %&% "__" %&% this_comp %&% "__de_upset.png",
                  width = 20 / 2.45, height = 10 / 2.45, units = "in", res = 300
                )
                print(p_upset)
                dev.off()
              }
            }
          }
        }
        # save
        message("Saving the DE results on table and on a RDS object")
        write.table(dres, outdir %&% this_ct %&% "__" %&% this_comp %&% "__deseq.txt",
          sep = "\t", row.names = FALSE
        )
        saveRDS(des, outdir %&% this_ct %&% "__" %&% this_comp %&% "__des.RDS")
      }
    }
  }
}
