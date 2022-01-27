#!/usr/bin/env Rscript
#################################################
## File name: cohort_aggregate.R
## Author: Matteo Carrara
## Date created: Aug 2021
## R Version: 4.1.0
##################################################

## TESTED ON THE FOLLOWING VERSIONS OF THE DEPENDENCIES:
## SingleCellExperiment 1.14.1
## reshape2 1.4.4
## optparse 1.6.6

## GENERAL:
## This script is the first part of a 3-part processing of single-cell data. Before moving to the analysis of the data, results from the single cell pipeline
## must be aggregated, so that instead of being separated by samples, it is a single composed element for a cohort of interest

### ADDITIONAL IMPORTANT INFORMATION:
## This script takes as input the results from the single cell pipeline. It assumes all input comes in the RDS format, you have one file per sample and you
## are processing each cohort separately. At the moment there is no way to process multiple cohort with a single run
## All options are required and must have a value set. In the future there might be a migration away from optparse and to a package allowing to
## automatically check if a required option has a value or not

## INPUT:
## --metadata A metadata file with two columns containing matches between RDS files and sampleID. Please refer to the section "DETAILS ON METADATA FILE" for more information.
## --groupvar The grouping feature used to aggregate the samples. It should be an existing column of the colData() of all objects
## --geneSelectionMethod The method used to select the genes. 'Union' is suggested and DEFAULT, because only a few genes are going to be shared between all samples
## --min_group_size The minimum number of cells in a cluster below which the cluster is filtered out. The default is 20
## --max_fractionMT The maximum allowed fraction of mitochondrial genes. The default is 0.5. If any input object does not have a colData() column called "fractionMT" it will be skipped with a warning
## - The directory in which the output must be stored. The script creates a subdirectory called sampleRDS
## - The level of verbosity of the script (0 = drop most messages, 1 = report activity regularly). DEFAULT = 1

## OUTPUT:
## - For each element of the grouping feature (e.g. for each celltype if grouping is done by celltype)
##   - <groupvar_value>_aggregate.<groupvar_name>.RDS contains the aggregated for the specific grouping variable
## - For each sample.
##   - sampleRDS/<sampleID>__<groupvar_name>.RDS is the count table for the specific sample
##   - samplesMeta/<sampleID>__<groupvar_name>_meta.txt is the sum of all counts in the sample for each value of "groupvar" within that same sample
## - aggregate.<groupvar_name>.txt is the aggregation of all sums available in folder "samplesMeta"

## REMARK on sample identifiers
## The sample identifiers are the only element that keeps track of the data throughout the analysis. It is very important to keep the identifiers
## consistent between input files. The two major elements that make use of the sample identifiers are the RDS files and the  clinical data table that can will be used in the next steps

## DETAILS ON METADATA FILE
# The metadata file stores the main information about the cohort files to aggregate
# The file must be a tab-delimited text table with 2 columns and no header
# The first column must contain the path to the RDS files to aggregate,, one file for each row
# The second column must contain the sampleID to use within the analysis for each specific file. Standard R variable rules apply (no numbers as first characters, aovid special characters and spaces)

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(reshape2)
  library(optparse)
})

# parse command line arguments
option_list <- list(
  make_option("--metadata", type = "character"),
  make_option("--groupvar", type = "character", help = "Feature to use to aggregate the samples. Must be an available field of colData of every input SCE object."),
  make_option("--geneSelectionMethod", type = "character", default = "union", help = "Either 'union' or 'intersect'. Decides if the script must work on the intersection or the union of features. Union is strongly recommended to increase the number of results."),
  make_option("--min_group_size", type = "integer", default = 20, help = "The minimum number of cells present in a group to be included in the analysis."),
  make_option("--max_fractionMT", type = "double", default = 0.5, help = "The maximum fraction of mitochondrial genes allowed in each cell. Must be between 0 and 1"),
  make_option("--outdir", type = "character", help = "Full path to output directory.")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

"%&%" <- function(a, b) paste(a, b, sep = "")
message(opt)
message("Setting up the input")
mymeta <- opt$metadata
mymeta <- read.delim(mymeta, header = FALSE, stringsAsFactors = FALSE)
rdsfn <- mymeta[, 1]
sampleID <- mymeta[, 2]

groupvar <- opt$groupvar
gene_select <- opt$geneSelectionMethod
min_group_size <- opt$min_group_size
max_fractionMT <- opt$max_fractionMT
outdir <- opt$outdir %&% "/"

number <- seq_along(rdsfn)

f_groupsum <- function(x) tapply(x, colData(sce)[, groupvar], sum, na.rm = TRUE)

message("Starting loading and grouping of the input files.\nThis process may take a while and use a lot of memory depending on the number of samples and their size.")

groups_all_meta <- list()
if (length(grep("samplesRDS", list.files(outdir))) == 0) {
  dir.create(outdir %&% "samplesRDS")
}
if (length(list.files(outdir %&% "samplesRDS")) != 0) {
  stop("The subdirectory 'samplesRDS' already exists in the output folder and is not empty")
}
if (length(grep("samplesMeta", list.files(outdir))) == 0) {
  dir.create(outdir %&% "samplesMeta")
}
if (length(list.files(outdir %&% "samplesMeta")) != 0) {
  stop("The subdirectory 'samplesMeta' already exists in the output folder and is not empty")
}

for (ii in number) {
  message(paste("Working on file", rdsfn[[ii]], "\tFile", ii, "of", length(number)))
  sce <- readRDS(rdsfn[ii])
  idx <- which(!sce$celltype_major %in% c("uncertain", "unknown") & sce$fractionMT < max_fractionMT)
  if (length(idx) == 0) {
    message("WARNING: sample " %&% rdsfn[ii] %&% "has no fractionMT column or no entry below the MT threshold. Skipping...")
    outfn <- outdir %&% "samplesRDS/" %&% sampleID[ii] %&% "__" %&% groupvar %&% ".RDS"
    m_group <- NA
    saveRDS(m_group, outfn)
    outfn <- outdir %&% "samplesMeta/" %&% sampleID[ii] %&% "__" %&% groupvar %&% "_meta.txt"
    write.table(NA, outfn, sep = "\t", row.names = FALSE)
    groups_all_meta[[ii]] <- NA
    next()
  }
  sce <- sce[, idx]
  if (!any(colnames(colData(sce)) %in% groupvar)) {
    stop(paste("Error: the value -", groupvar, "- for 'groupvar' is not available in the sce object", rdsfn[ii]))
  }
  colData(sce)[, groupvar] <- as.character(colData(sce)[, groupvar])
  groups_all <- melt(table(colData(sce)[, groupvar]),
    as.is = T,
    value.name = "counts", varnames = "cluster"
  )
  groups_all$type <- groupvar
  id_keep <- groups_all$cluster[groups_all$counts >= min_group_size]
  groups_all <- groups_all[groups_all$cluster %in% id_keep, , drop = FALSE]
  m_group <- apply(counts(sce), 1, f_groupsum)

  # We need to add [1] here or it throws a warning
  if (class(m_group)[1] == "numeric") {
    m_group <- matrix(m_group, nrow = 1)
    rownames(m_group) <- unique(colData(sce)[, groupvar])[1]
    colnames(m_group) <- rownames(counts(sce))
  }

  outfn <- outdir %&% "samplesRDS/" %&% sampleID[ii] %&% "__" %&% groupvar %&% ".RDS"
  m_group <- t(m_group[id_keep, , drop = FALSE])
  saveRDS(m_group, outfn)
  outfn <- outdir %&% "samplesMeta/" %&% sampleID[ii] %&% "__" %&% groupvar %&% "_meta.txt"
  write.table(groups_all, outfn, sep = "\t", row.names = FALSE)
  groups_all_meta[[ii]] <- groups_all
}

group_list <- list()
meta_df <- data.frame(NULL)

for (ii in seq_along(rdsfn)) {
  if (is.na(unlist(groups_all_meta[[ii]])[1])) {
    next()
  }
  this_meta <- groups_all_meta[[ii]]
  this_meta <- this_meta[this_meta$type == groupvar, , drop = FALSE]
  this_meta$sample <- sampleID[ii]
  meta_df <- rbind(meta_df, this_meta)
  this_groups <- this_meta$cluster
  m_group <- readRDS(outdir %&% "samplesRDS/" %&% sampleID[ii] %&% "__" %&% groupvar %&% ".RDS")
  for (jj in seq_along(this_groups)) {
    x <- data.frame(m_group[, jj, drop = FALSE])
    names(x) <- this_meta$sample[jj]
    rownames(x)[1:5]
    if (this_groups[jj] %in% names(group_list)) {
      if (gene_select == "intersect") {
        genes_keep <- intersect(
          rownames(group_list[[this_groups[jj]]]),
          rownames(x)
        )
        group_list[[this_groups[jj]]] <- cbind(
          group_list[[this_groups[jj]]][genes_keep, , drop = FALSE],
          x[genes_keep, , drop = FALSE]
        )
      } else if (gene_select == "union") {
        old_genes <- rownames(group_list[[this_groups[jj]]])
        all_genes <- union(old_genes, rownames(x))
        new_genes <- setdiff(all_genes, old_genes)
        m_empty <- matrix(0, nrow = length(new_genes), ncol = ncol(group_list[[this_groups[jj]]]))
        colnames(m_empty) <- colnames(group_list[[this_groups[jj]]])
        rownames(m_empty) <- new_genes
        group_list[[this_groups[jj]]] <- rbind(group_list[[this_groups[jj]]], m_empty)
        m_new <- matrix(0, nrow = length(all_genes), ncol = 1)
        rownames(m_new) <- all_genes
        colnames(m_new) <- names(x)
        m_new[rownames(x), ] <- x[, 1]
        group_list[[this_groups[jj]]] <- cbind(group_list[[this_groups[jj]]], m_new)
      }
    } else {
      lsize <- length(group_list)
      group_list <- c(group_list, list(x))
      names(group_list)[lsize + 1] <- this_groups[jj]
    }
  }
}

write.table(meta_df[
  order(meta_df$cluster, meta_df$sample),
  c("cluster", "sample", "counts")
],
file = outdir %&% "aggregate." %&% groupvar %&% ".txt",
sep = "\t", row.names = FALSE
)
for (ii in seq_len(length(group_list))) {
  fn <- outdir %&% names(group_list)[ii] %&% "_aggregate." %&% groupvar %&% ".RDS"
  saveRDS(group_list[[ii]], fn)
}
