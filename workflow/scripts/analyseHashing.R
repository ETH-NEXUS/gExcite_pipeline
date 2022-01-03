#############
## Script to analyse hashed scData.
## Produces barcode list / tag that can be used in the follow up analysis.
#############

# Input: Script takes the results of a Cellranger Run and a Citeseq Run as an argument.
# Analysis: The cells (barcodes) from both runs are intersected and analysed using Seurat.
# Output: One file containing barcodes per tag used in the citeseq run + Overview / QC Plots

# Linda Grob


library(optparse)
library(Seurat)
library(cowplot)
library(ggplot2)
library(viridisLite)
library(parallelDist)
library(scales)
library(uwot)
library(patchwork)
library(plyr)
library(glue)

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")


option_list <- list(
  make_option("--citeseq_in", type = "character", help = "Path to the input folder that contains the citeseq results."),
  make_option("--adt_barcodes_in", type = "character", help = "Path to the input files that contains the adt barcodes."),
  make_option("--quantile_threshold", type = "character", help = "Quantile threshold for HTODemux() function. E.g. 0.99"),
  make_option("--sampleTagMap", type = "character", help = "Map associating barcodes, tags and Samplenames."),
  make_option("--normalisation",
    type = "character",
    help = "Parameter for normalisation. Can be '1' (normalisation across features), '2' (normalisation across cells), or '1,2' (both normalisations combined)."
  ),
  make_option("--normalisation_downstream",
    type = "character",
    help = "If normalisation = '1,2' the normalised counts of only one option are kept. This option can be '1' or '2'."
  ),
  make_option("--save_negatives",
    type = "logical", default = TRUE,
    help = "If normalisation '1,2' is chosen this parameter decides if a given tag stays for the cell even if using the second normalisation method a 'Negative' is predicted. Can be 'FALSE' or 'TRUE', default is TRUE."
  ),
  make_option("--output_prefix", type = "character", help = "Path that is used as prefix for files.")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# check option --save_negatives
print("Option '--save_negatives' is:")
print(opt$save_negatives)

# create output directory
print("Creating main output dir:")
print(dirname(opt$output_prefix))
dir.create(dirname(opt$output_prefix))

# get output prefix for files
prefix_only <- gsub(pattern = "^.*/([^/]+)$", "\\1", opt$output_prefix)
# create subdirectory for plots
print("Creating subdirectory for plots:")
outdir_plots <- paste0(dirname(opt$output_prefix), "/plots/")
print(outdir_plots)
dir.create(outdir_plots)
plot_prefix <- paste0(outdir_plots, prefix_only)
print(plot_prefix)

# Read in SampleMap
sampleTagMap <- read.csv(opt$sampleTagMap, sep = ",", header = FALSE)
colnames(sampleTagMap) <- c("barcode", "tagName", "sample")
print(sampleTagMap)

# Read in parameter
quantile_threshold <- as.numeric(opt$quantile_threshold)
print(paste0("Quantile threshold: ", quantile_threshold))

# Read in data
CiteSeq <- Read10X(data.dir = opt$citeseq_in, gene.column = 1)
CS.HTO <- CreateAssayObject(counts = CiteSeq)
CS.HTO <- as.matrix(CS.HTO[-nrow(CS.HTO), ])
glue(str(CS.HTO))
glue("Cell barcodes of cite-seq-count CS.HTO")
print(head(colnames(CS.HTO)))

# The newest version of CiteSeq is able to run multiple lanes at the same time. --> To do adapt snake make pipeline and this script.

CellRanger_results.ADT <- Read10X(data.dir = opt$adt_barcodes_in)
glue(str(CellRanger_results.ADT))
glue("Cell barcodes of CellRanger_results.ADT")
print(head(colnames(CellRanger_results.ADT)))

table(grepl("-1", colnames(CellRanger_results.ADT)))
cat("\n\n\nStop if not all cellranger cell barcodes have a '-1' as suffix\n\n")
stopifnot(length(colnames(CellRanger_results.ADT)) == table(grepl("-1", colnames(CellRanger_results.ADT))))
colnames(CellRanger_results.ADT) <- gsub("-1", "", colnames(CellRanger_results.ADT))


# Only continue with barcodes found in citeseq and cellranger analysis
joint.bcs <- intersect(colnames(CellRanger_results.ADT), colnames(CS.HTO))
print(dim(CS.HTO))
CS.HTO <- as.matrix(CS.HTO[, joint.bcs])
rownames(CS.HTO) <- gsub("[TCGA-]", "", rownames(CS.HTO))

CellRanger_results.ADT <- CellRanger_results.ADT[, joint.bcs]


#################################################################################
### Perform count normalisation either across cells, across features, or both ###
#################################################################################

# read in normalisation parameters
norm_methods <- strsplit(opt$normalisation, split = ",")
norm_methods <- unlist(norm_methods)
norm_methods <- as.integer(norm_methods)
str(norm_methods)
norm_downstream <- as.integer(opt$normalisation_downstream)
str(norm_downstream)
stopifnot(norm_downstream %in% c(1, 2))

# print out chosen normalisation option
if (all.equal(norm_methods, 1) == TRUE) {
  print("Normalisation method: Across features")
} else if (all.equal(norm_methods, 2) == TRUE) {
  print("Normalisation method: Across cells")
} else if (all.equal(norm_methods, c(1, 2)) == TRUE) {
  print("Normalisation method: Combine normalisation across cells and across features")
} else {
  stop("The parameter 'normalisation' was not set correctly.")
}

# initiate list for seurat objects
results_demux <- list()
# perform HTO demux
for (method in norm_methods) {
  print(paste0("Normalising with margin = ", method))
  # Store data in Seurat object
  my_sample <- CreateSeuratObject(counts = CellRanger_results.ADT)
  print(my_sample)
  my_sample <- NormalizeData(my_sample)
  my_sample <- FindVariableFeatures(my_sample)
  my_sample <- ScaleData(my_sample)
  print(str(CS.HTO))
  my_sample[["HTO"]] <- CreateAssayObject(counts = CS.HTO)
  my_sample <- NormalizeData(my_sample, assay = "HTO", normalization.method = "CLR", margin = method)
  my_sample <- ScaleData(my_sample, assay = "HTO")
  ## Demux Sample
  ## The HTODemux function seems to produce an error if either to many cells are found in the sample or the data is too noisy.
  ## Alternative function is used in that case - however then only one list containing all tags is created together with an error message file.
  ## If it becomes clear why this happens, potentially adapt the code here.
  function_error <- FALSE

  tryCatch(
    {
      my_sample <- HTODemux(my_sample, assay = "HTO", positive.quantile = quantile_threshold)
    },
    warning = function(w) {
      writeLines("The HTODemux function produced a warning.", paste(opt$output_prefix, "warning.txt"), sep = ".")
    },
    error = function(e) {
      function_error <<- TRUE
      writeLines("The HTODemux function produced an error. Used MULTIseqDemux instead.", paste(trimws(opt$output_prefix), "FUNCTION_ERROR.txt", sep = "."))
      my_sample <<- MULTIseqDemux(my_sample, assay = "HTO", autoThresh = TRUE)
    }
  )

  # Generate Count Table based on the maxID Tag that Seurat assigns
  ## "Raw Count" shows number of cells / "best-hit-hashtags" without taking singlets / doublets into account.
  if (function_error == FALSE) {
    print("Hashing analysis with HTODemux was successful.")
    print("Proceeding with plotting and writing results to files.")
  } else {
    # In case HTO demux fails only one list containing all tags is produced at the moment. Either remove that part or add addditional ouput lists.
    my_sample.subset <- subset(my_sample, idents = "Negative", invert = TRUE)
    write.table(colnames(my_sample.subset[which(my_sample.subset$MULTI_ID != "Doublet")]), paste(opt$output_prefix, "barcodes_singlets.txt", sep = "."), row.names = FALSE, col.names = FALSE, quote = FALSE)
  }

  # Proceeding if HTODemux did not fail
  table(my_sample$HTO_maxID)
  write.table(table(my_sample$HTO_maxID),
    file = paste(opt$output_prefix, method, "normalisation_raw_tag_counts.txt", sep = "."),
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
  )

  ## Cell Counts / hashtag with Doublets & Negative cells treated special.
  table(my_sample$hash.ID)
  write.table(table(my_sample$hash.ID),
    file = paste(opt$output_prefix, method, "normalisation_final_tag_counts.txt", sep = "."),
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
  )

  # RIDGEPLOT
  # Group Cells based on the highest hashtag signal.
  print("Plotting RidgePlot")
  RidgePlot(my_sample, assay = "HTO", features = rownames(my_sample[["HTO"]]), ncol = 2)
  rigdePlot <- RidgePlot(my_sample, assay = "HTO", features = rownames(my_sample[["HTO"]]), ncol = 2)
  ggplot2::ggsave(filename = paste(opt$output_prefix, method, "normalisation_ridgeplot.png", sep = "."), plot = rigdePlot)

  # HTOHeatmap
  print("Plotting Heatmap")
  hto_heatmap <- HTOHeatmap(my_sample, assay = "HTO", raster = FALSE)
  ggplot2::ggsave(filename = paste(opt$output_prefix, method, "normalisation_htoheatmap.png", sep = "."), plot = hto_heatmap)

  # put resulting Seurat object into list
  method_string <- paste0("method_", method)
  results_demux[method_string] <- my_sample
}

print("List of HTO_Demux results:")
print(str(results_demux, max.level = 2))

#####
# if more than one HTO demux object, then compare and filter
if (length(results_demux) == 1) {
  print("Only one normalisation option was chosen. Continuing with the plots.")
} else if (length(results_demux) == 2) {
  print("Both normalisation options were chosen. Continuing with filtering the object.")
  cellIDtable_1 <- results_demux[[1]][[]]
  stopifnot(all.equal(results_demux[["method_1"]], results_demux[[1]]))
  cellIDtable_1$barcodes <- rownames(cellIDtable_1)
  cellIDtable_2 <- results_demux[[2]][[]]
  stopifnot(all.equal(results_demux[["method_2"]], results_demux[[2]]))
  cellIDtable_2$barcodes <- rownames(cellIDtable_2)
  both_results <- dplyr::inner_join(cellIDtable_1, cellIDtable_2, by = "barcodes", suffix = c("_1", "_2"))

  # All columns that are the same for both normalisations are reduced to one column
  print("")
  print("Reduce identical columns:")
  for (dataCol in c("orig.ident", "nCount_RNA", "nFeature_RNA", "nCount_HTO", "nFeature_HTO")) {
    print(dataCol)
    col1 <- paste0(dataCol, "_1")
    col2 <- paste0(dataCol, "_2")
    stopifnot(all.equal(both_results[[col1]], both_results[[col2]]))
    both_results[dataCol] <- both_results[col1]
    both_results[col1] <- NULL
    both_results[col2] <- NULL
  }
  # have HTO_classification columns as characters
  both_results$HTO_classification_1 <- as.character(both_results$HTO_classification_1)
  both_results$HTO_classification_2 <- as.character(both_results$HTO_classification_2)

  # if a 'Negative' and a doublet were detected the tag from the normalisation chosen for downstream analyses will be chosen.
  classification_downstream <- paste0("HTO_classification_", norm_downstream)

  #####     Compare results in columns 'HTO_classification_1' and 'HTO_classification_2'     #####
  HTO_class_1 <- both_results$HTO_classification_1
  HTO_class_2 <- both_results$HTO_classification_2
  HTO_class <- rep("Error", length(HTO_class_1))

  # initiate counts for the different combinations
  count_identical <- 0
  count_neg_doub <- 0
  count_neg <- 0
  count_neg_saved <- 0
  count_doub <- 0
  count_sing_sing <- 0
  count_error <- 0

  for (i in seq_along(HTO_class_1)) {
    ### tags identical
    if (HTO_class_1[i] == HTO_class_2[i]) {
      HTO_class[i] <- HTO_class_1[i]
      count_identical <- count_identical + 1
    }
    ### 'Negative' detected
    else if (HTO_class_1[i] == "Negative" | HTO_class_2[i] == "Negative") {
      # with doublet and Negative take the result of the normalisation chosen for downstream analyses
      if (grepl("_", HTO_class_1[i]) | grepl("_", HTO_class_2[i])) {
        HTO_class[i] <- both_results[[classification_downstream]][i]
        count_neg_doub <- count_neg_doub + 1
      }
      # save Negatives is false
      else if (!opt$save_negatives) {
        HTO_class[i] <- "Negative"
        count_neg <- count_neg + 1
      }
      # save Negatives is true and 'HTO_classification_1[i]' has singlet tag
      else if (opt$save_negatives & grepl("^tag[^_]+$", HTO_class_1[i])) {
        HTO_class[i] <- HTO_class_1[i]
        count_neg_saved <- count_neg_saved + 1
      }
      # save Negatives is true and 'HTO_classification_2[i]' has singlet tag
      else if (opt$save_negatives & grepl("^tag[^_]+$", HTO_class_2[i])) {
        HTO_class[i] <- HTO_class_2[i]
        count_neg_saved <- count_neg_saved + 1
      }
    }
    ### doublet detected and other result is either different doublet or singlet
    else if (grepl("_", HTO_class_1[i]) | grepl("_", HTO_class_2[i])) {
      HTO_class[i] <- "Doublet"
      count_doub <- count_doub + 1
    }
    ### two different singlets detected
    else if (grepl("^tag[^_]+$", HTO_class_1[i]) & grepl("^tag[^_]+$", HTO_class_2[i])) {
      HTO_class[i] <- "Doublet"
      count_sing_sing <- count_sing_sing + 1
    }
    ### if any case is not covered by the above insert "Error"
    else {
      HTO_class[i] <- "Error"
      count_error <- count_error + 1
    }
  }
  stopifnot(sum(HTO_class == "Error") == 0)

  # print out results
  cat("\n")
  cat("Results of comparing both normalisation approaches:\n")
  cat("Results identical: ", count_identical, "\n")
  cat("Negative and doublet: ", count_neg_doub, "\n")
  cat("Negatives not saved: ", count_neg, "\n")
  cat("Negatives saved: ", count_neg_saved, "\n")
  cat("Doublets (either two different doublets or doublet and singlet): ", count_doub, "\n")
  cat("Two different Singlets: ", count_sing_sing, "\n")
  cat("Errors: ", count_error, "\n\n")

  # add resulting HTO_class vector as "HTO_classification" to table
  both_results$HTO_classification <- HTO_class

  # infer column "HTO_classification" from "HTO_classification"
  both_results$HTO_classification.global <- ifelse(both_results$HTO_classification == "Negative",
    # if Negative in HTO_classification then also in global
    "Negative",
    ifelse(both_results$HTO_classification == "Doublet",
      # if Doublet in HTO_classification then also in global
      "Doublet",
      ifelse(grepl("_", both_results$HTO_classification),
        # if Doublet in HTO_classification then also in global
        "Doublet",
        ifelse(grepl("^tag[^_]+$", both_results$HTO_classification),
          # if any of the tags is in HTO_classification then it is Singlet in global
          "Singlet",
          # Give error if anything is not covered by the cases above
          "Error"
        )
      )
    )
  )
  stopifnot(sum(both_results$HTO_classification.global == "Error") == 0)

  # infer column "hash.ID" from "HTO_classification"
  both_results$hash.ID <- ifelse(both_results$HTO_classification == "Negative",
    # if Negative in HTO_classification then also in hash.ID
    "Negative",
    ifelse(both_results$HTO_classification == "Doublet",
      # if Doublet in HTO_classification then also in hash.ID
      "Doublet",
      ifelse(grepl("_", both_results$HTO_classification),
        # if Doublet in HTO_classification then also in hash.ID
        "Doublet",
        ifelse(grepl("^tag[^_]+$", both_results$HTO_classification),
          # if any of the tags is in HTO_classification then it is the same in hash.ID)
          both_results$HTO_classification,
          # Give error if anything is not covered by the cases above
          "Error"
        )
      )
    )
  )
  stopifnot(sum(both_results$hash.ID == "Error") == 0)

  # write full table into file
  filename_both_results <- paste0(opt$output_prefix, ".both_normalisation_results_HTODemux.tsv")
  print(filename_both_results)
  write.table(
    x = both_results, file = filename_both_results,
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
  )

  # into the column "HTO_margin" write the results of the normalisation chosen for downstream analyses
  HTO_margin_downstream <- paste0("HTO_margin_", norm_downstream)
  both_results$HTO_margin <- both_results[[HTO_margin_downstream]]

  # have one of the two seurat objects as "my_sample"
  # which of the normalised counts are in the object should not matter
  # The UMAP is calculated on the raw counts and the normalised counts are not saved
  # in case the counts are needed for any downstream analyses a parameter can be set to choose one of the two normalisation options
  my_sample <- results_demux[[norm_downstream]]

  # in my_sample remove all "old" columns to not confuse them with the newly inferred
  all_old_columns <- names(my_sample[[]])
  cat("\n")
  print("Remove old columns:")
  for (old_col in all_old_columns) {
    print(old_col)
    my_sample[[old_col]] <- NULL
  }
  # insert the new hash tag data columns that were inferred from both normalisations
  cat("\n")
  print("Insert columns inferred from both normalisations:")
  for (new_col in names(both_results)) {
    print(new_col)
    my_sample[[new_col]] <- both_results[[new_col]]
  }
  stopifnot(rownames(my_sample[[]]) == my_sample$barcodes)
  # make sure the data format is the same as originally
  my_sample$HTO_classification <- as.factor(my_sample$HTO_classification)
  my_sample$HTO_classification.global <- as.factor(my_sample$HTO_classification.global)
  my_sample$hash.ID <- as.factor(my_sample$hash.ID)
  my_sample$HTO_margin <- as.numeric(my_sample$HTO_margin)

  # replace two meta data columns of the downstream option with the default names
  # they are required later in this script
  maxID <- paste0("HTO_maxID_", norm_downstream)
  secondID <- paste0("HTO_secondID_", norm_downstream)
  newNames <- gsub(maxID, "HTO_maxID", names(my_sample[[]]))
  newNames <- gsub(secondID, "HTO_secondID", newNames)
  names(my_sample@meta.data) <- newNames
} else {
  stop("Something went wrong with the count normalisation. Please check the input.")
}


## Overview of number of counts per tag, including Doublets and Negatives
## If both normalisation were applied, this shows the final results
table(my_sample$hash.ID)
filename_final_tag_counts <- paste(opt$output_prefix, "final_tag_counts.txt", sep = ".")
write.table(table(my_sample$hash.ID), file = filename_final_tag_counts, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


cat("\n")
print("Start plotting.")
print("Seurat object:")
print(my_sample)

################################################################
###   Format data for plotting histograms and scatterplots   ###
################################################################

# format sparse HTO raw count matrix to data frame
hto_counts <- GetAssayData(object = my_sample, assay = "HTO", slot = "counts")
hto_counts <- as.matrix(hto_counts)
hto_counts <- t(hto_counts)
hto_counts <- as.data.frame(hto_counts)
hto_counts$barcodes <- rownames(hto_counts)
# get cell metadata into data frame
cell_data <- as.data.frame(my_sample[[]])
cell_data$barcodes <- rownames(cell_data)
stopifnot(rownames(hto_counts) == rownames(cell_data))
# join cell meta data and hto counts
cellData_htoCounts <- dplyr::inner_join(cell_data, hto_counts, by = "barcodes")
all_hash.IDs <- sort(unique(cellData_htoCounts$hash.ID))

# Infer minimum raw count per tag
all_tags <- all_hash.IDs[!all_hash.IDs %in% c("Doublet", "Negative")]
all_tags <- as.character(all_tags)
# get minimum value
inferred_minCounts <- numeric()
for (tag in all_tags) {
  tag_cellData_htoCounts <- subset(cellData_htoCounts, hash.ID == tag)
  stopifnot(length(unique(tag_cellData_htoCounts$hash.ID)) == 1)
  tag_minCount <- min(tag_cellData_htoCounts[[tag]])
  print(paste0("minCount of tag ", tag, ": ", tag_minCount))
  inferred_minCounts[tag] <- tag_minCount
}
# Write inferred minCount into table
minCounts_out <- as.data.frame(inferred_minCounts)
names(minCounts_out) <- "minimum_raw_count"
minCounts_out$hashtag <- names(inferred_minCounts)
filename_minCounts <- paste(opt$output_prefix, "minimum_raw_count_per_tag.tsv", sep = ".")
print(filename_minCounts)
write.table(
  x = minCounts_out, file = filename_minCounts,
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
)


###########################
###   Plot histograms   ###
###########################

# have fixed colours for the histogram fill parameter
# E.g. tag5 should have always the same colour
colours_hash.IDs <- hue_pal()(length(all_hash.IDs))
names(colours_hash.IDs) <- all_hash.IDs

# colour setting for histograms based on cellData_htoCounts
cellData_htoCounts_colours <- match(cellData_htoCounts$hash.ID, names(colours_hash.IDs))

# histogram of HTO_margin
histo_margin <- ggplot(cellData_htoCounts, aes(x = HTO_margin, fill = hash.ID)) +
  geom_histogram(binwidth = 0.1) +
  scale_fill_manual(name = "hash.ID", values = colours_hash.IDs[cellData_htoCounts_colours]) +
  ggtitle("HTO_margin (difference between signals for hash.maxID and hash.secondID)")
histo_margin
file_name <- paste(plot_prefix, "histo.hto_margin.png", sep = ".")
print(file_name)
ggplot2::ggsave(filename = file_name, plot = histo_margin, width = 22, height = 17, units = "cm")

# histogram of raw counts
histo_htoCounts <- ggplot(cellData_htoCounts, aes(x = nCount_HTO, fill = hash.ID)) +
  geom_histogram(binwidth = 30) +
  coord_cartesian(xlim = c(0, quantile(cellData_htoCounts$nCount_HTO, .99))) +
  scale_fill_manual(name = "hash.ID", values = colours_hash.IDs[cellData_htoCounts_colours]) +
  ggtitle("Raw counts")
histo_htoCounts
file_name <- paste(plot_prefix, "histo.hto_counts.png", sep = ".")
print(file_name)
ggplot2::ggsave(filename = file_name, plot = histo_htoCounts, width = 19, height = 17, units = "cm")


###   Plot histograms only Singlets   ###
data_only_singlets <- subset(cellData_htoCounts, HTO_classification.global == "Singlet")
# colour setting for histograms based on data_only_singlets
singlets_colours <- match(data_only_singlets$hash.ID, names(colours_hash.IDs))

# histogram of margins
singlet_histo_margin <- ggplot(data_only_singlets, aes(x = HTO_margin, fill = hash.ID)) +
  geom_histogram(binwidth = 0.1) +
  scale_fill_manual(name = "hash.ID", values = colours_hash.IDs[singlets_colours]) +
  ggtitle("HTO_margin (difference between signals for hash.maxID and hash.secondID)")
singlet_histo_margin
file_name <- paste(plot_prefix, "singlet_histo.hto_margin.png", sep = ".")
print(file_name)
ggplot2::ggsave(filename = file_name, plot = singlet_histo_margin, width = 22, height = 17, units = "cm")

# histogram of raw counts
singlet_histo_htoCounts <- ggplot(data_only_singlets, aes(x = nCount_HTO, fill = hash.ID)) +
  geom_histogram(binwidth = 30) +
  coord_cartesian(xlim = c(0, quantile(data_only_singlets$nCount_HTO, .99))) +
  scale_fill_manual(name = "hash.ID", values = colours_hash.IDs[singlets_colours]) +
  ggtitle("Raw counts")
singlet_histo_htoCounts
file_name <- paste(plot_prefix, "singlet_histo.hto_counts.png", sep = ".")
print(file_name)
ggplot2::ggsave(filename = file_name, plot = singlet_histo_htoCounts, width = 19, height = 17, units = "cm")

# give out one plot with all histograms
all_histo_plots <- list(
  histo_margin, histo_htoCounts,
  singlet_histo_margin, singlet_histo_htoCounts
)
my_grid <- cowplot::plot_grid(plotlist = all_histo_plots, ncol = 2)
file_all_histo <- paste(opt$output_prefix, "all_histograms.png", sep = ".")
print(file_all_histo)
ggplot2::ggsave(filename = file_all_histo, plot = my_grid, width = 40, height = 15 * ceiling(length(all_histo_plots) / 2), units = "cm")


#############################
###   Plot Scatterplots   ###
#############################
cat("plotting scatterplots\n")

# get all combinations of two hashtags
all_hashtags <- rownames(my_sample[["HTO"]])
pairs_tags <- combn(all_hashtags, 2, simplify = FALSE)
# give out scatterplots for each pair of hash tag oligos (hto)
for (pair in pairs_tags) {
  print("Pair of hash tag oligos:")
  print(pair)
  tag_y <- pair[1]
  tag_x <- pair[2]
  pair_data <- subset(cellData_htoCounts, HTO_maxID == tag_y | HTO_maxID == tag_x)
  stopifnot(sort(as.character(unique(pair_data$HTO_maxID))) == sort(pair))
  # rename Singlet to include the identity of the singlet
  hashing_results <- ifelse(pair_data$HTO_classification.global == "Doublet", "Doublet",
    ifelse(pair_data$HTO_classification.global == "Negative", "Negative",
      ifelse(pair_data$hash.ID == tag_x, paste0(tag_x, "_Singlet"),
        paste0(tag_y, "_Singlet")
      )
    )
  )
  print(str(hashing_results))
  print(table(hashing_results))
  stopifnot(length(table(hashing_results)) == 4)
  # find for each pair the maximum for x and y axis
  # y and x have the same axis limit that is the mean of the 0.999 quantile of both
  limit_x <- quantile(pair_data[[tag_x]], .99)
  limit_y <- quantile(pair_data[[tag_y]], .99)
  axis_max <- mean(c(limit_x, limit_y))

  ########################
  ####    PLOTTING    ####
  ########################
  # scatterplot with all cells showing the counts of tag_x and tag_y
  # colour is detailed HTO_classification
  scatter_plot_all <- ggplot() +
    geom_point(shape = 1, stroke = .9, aes(
      x = cellData_htoCounts[[tag_x]],
      y = cellData_htoCounts[[tag_y]],
      color = cellData_htoCounts$HTO_classification
    )) +
    coord_cartesian(
      xlim = c(0, 500),
      ylim = c(0, 500)
    ) +
    xlab(names(pair_data[tag_x])) +
    ylab(names(pair_data[tag_y])) +
    scale_color_discrete(name = "HTO_classification") +
    geom_vline(xintercept = inferred_minCounts[tag_x], color = "blue3") +
    annotate("text",
      x = inferred_minCounts[tag_x], y = 0,
      label = paste0(names(inferred_minCounts[tag_x]), ": ", inferred_minCounts[tag_x]),
      hjust = -0.3, vjust = 2
    ) +
    geom_hline(yintercept = inferred_minCounts[tag_y], color = "blue3") +
    annotate("text",
      x = 0, y = inferred_minCounts[tag_y],
      label = paste0(names(inferred_minCounts[tag_y]), ": ", inferred_minCounts[tag_y]),
      vjust = -0.9, hjust = 0.2
    ) +
    ggtitle(paste0("All cells, raw counts of ", tag_x, " and ", tag_y))
  scatter_plot_all
  file_name <- paste(plot_prefix, "scatterplot_raw_counts.all_cells", tag_x, tag_y, "png", sep = ".")
  print(file_name)
  ggplot2::ggsave(filename = file_name, plot = scatter_plot_all, width = 19, height = 17, units = "cm")


  # scatterplot with only cells with the highest signal in tag_x or tag_y
  # colour is detailed HTO_classification
  scatter_plot_doublet_details <- ggplot() +
    geom_point(shape = 1, stroke = .9, aes(
      x = pair_data[[tag_x]],
      y = pair_data[[tag_y]],
      color = pair_data$HTO_classification
    )) +
    coord_cartesian(
      xlim = c(0, limit_x),
      ylim = c(0, limit_y)
    ) +
    xlab(names(pair_data[tag_x])) +
    ylab(names(pair_data[tag_y])) +
    scale_color_discrete(name = "HTO_classification") +
    geom_vline(xintercept = inferred_minCounts[tag_x], color = "blue3") +
    annotate("text",
      x = inferred_minCounts[tag_x], y = 0,
      label = paste0(names(inferred_minCounts[tag_x]), ": ", inferred_minCounts[tag_x]),
      hjust = -0.3, vjust = 2
    ) +
    geom_hline(yintercept = inferred_minCounts[tag_y], color = "blue3") +
    annotate("text",
      x = 0, y = inferred_minCounts[tag_y],
      label = paste0(names(inferred_minCounts[tag_y]), ": ", inferred_minCounts[tag_y]),
      vjust = -0.9, hjust = 0.2
    ) +
    ggtitle(paste0("Only cells with HTO_maxID == ", tag_x, " | ", tag_y))
  scatter_plot_doublet_details
  file_name <- paste(plot_prefix, "scatterplot_raw_counts.doublet_details", tag_x, tag_y, "png", sep = ".")
  print(file_name)
  ggplot2::ggsave(filename = file_name, plot = scatter_plot_doublet_details, width = 19, height = 17, units = "cm")


  # scatterplot with only cells with the highest signal in tag_x or tag_y
  # both axis limited to the 0.999 quantile, respectively
  # colouring is Doublet, Negative, or tag_x or tag_y
  scatter_plot <- ggplot() +
    geom_point(shape = 1, stroke = .9, aes(
      x = pair_data[[tag_x]],
      y = pair_data[[tag_y]],
      color = hashing_results
    )) +
    coord_cartesian(
      xlim = c(0, limit_x),
      ylim = c(0, limit_y)
    ) +
    scale_colour_brewer(palette = "Dark2") +
    xlab(names(pair_data[tag_x])) +
    ylab(names(pair_data[tag_y])) +
    geom_vline(xintercept = inferred_minCounts[tag_x], color = "blue3") +
    annotate("text",
      x = inferred_minCounts[tag_x], y = 0,
      label = paste0(names(inferred_minCounts[tag_x]), ": ", inferred_minCounts[tag_x]),
      hjust = -0.3, vjust = 2
    ) +
    geom_hline(yintercept = inferred_minCounts[tag_y], color = "blue3") +
    annotate("text",
      x = 0, y = inferred_minCounts[tag_y],
      label = paste0(names(inferred_minCounts[tag_y]), ": ", inferred_minCounts[tag_y]),
      vjust = -0.9, hjust = 0.2
    ) +
    ggtitle(paste0("Only cells with HTO_maxID == ", tag_x, " | ", tag_y))
  scatter_plot
  file_name <- paste(plot_prefix, "scatterplot_raw_counts", tag_x, tag_y, "png", sep = ".")
  print(file_name)
  ggplot2::ggsave(filename = file_name, plot = scatter_plot, width = 19, height = 17, units = "cm")


  # scatterplot with only cells with the highest signal in tag_x or tag_y
  # both axis have same length, which is mean of both 0.999 quantiles
  # colouring is tag_x highlighted
  scatter_plot2 <- ggplot() +
    geom_point(shape = 1, stroke = .9, aes(
      x = pair_data[[tag_x]],
      y = pair_data[[tag_y]],
      color = hashing_results
    )) +
    coord_cartesian(
      xlim = c(0, axis_max),
      ylim = c(0, axis_max)
    ) +
    scale_colour_manual(values = c("grey", "grey", "grey", "red2")) +
    xlab(names(pair_data[tag_x])) +
    ylab(names(pair_data[tag_y])) +
    # to have all cells of interest on the top
    geom_point(
      shape = 1, stroke = .9, aes(
        x = pair_data[[tag_x]],
        y = pair_data[[tag_y]],
        color = hashing_results
      ),
      color = ifelse(grepl(tag_x, hashing_results), "red2", NA)
    ) +
    geom_vline(xintercept = inferred_minCounts[tag_x], color = "darkorchid4") +
    annotate("text",
      x = inferred_minCounts[tag_x], y = 0,
      label = paste0(names(inferred_minCounts[tag_x]), ": ", inferred_minCounts[tag_x]),
      hjust = -0.3, vjust = 2, color = "darkorchid4"
    ) +
    geom_hline(yintercept = inferred_minCounts[tag_y], color = "blue4") +
    annotate("text",
      x = 0, y = inferred_minCounts[tag_y],
      label = paste0(names(inferred_minCounts[tag_y]), ": ", inferred_minCounts[tag_y]),
      vjust = -0.9, hjust = 0.2, color = "blue4"
    ) +
    ggtitle(paste0("Only cells with HTO_maxID == ", tag_x, " or HTO_maxID == ", tag_y))
  scatter_plot2
  file_name2 <- paste(plot_prefix, "scatter_raw_counts", paste0(tag_x, "-highlighted"), tag_x, tag_y, "png", sep = ".")
  print(file_name2)
  ggplot2::ggsave(filename = file_name2, plot = scatter_plot2, width = 19, height = 17, units = "cm")


  # scatterplot with only cells with the highest signal in tag_x or tag_y
  # both axis have same length, which is mean of both 0.999 quantiles
  # colouring is tag_y highlighted
  scatter_plot3 <- ggplot() +
    geom_point(shape = 1, stroke = .9, aes(
      x = pair_data[[tag_x]],
      y = pair_data[[tag_y]],
      color = hashing_results
    )) +
    coord_cartesian(
      xlim = c(0, axis_max),
      ylim = c(0, axis_max)
    ) +
    scale_colour_manual(values = c("grey", "grey", "red2", "grey")) +
    xlab(names(pair_data[tag_x])) +
    ylab(names(pair_data[tag_y])) +
    # to have all cells of interest on the top
    geom_point(
      shape = 1, stroke = .9, aes(
        x = pair_data[[tag_x]],
        y = pair_data[[tag_y]],
        color = hashing_results
      ),
      color = ifelse(grepl(tag_y, hashing_results), "red2", NA)
    ) +
    geom_vline(xintercept = inferred_minCounts[tag_x], color = "darkorchid4") +
    annotate("text",
      x = inferred_minCounts[tag_x], y = 0,
      label = paste0(names(inferred_minCounts[tag_x]), ": ", inferred_minCounts[tag_x]),
      hjust = -0.3, vjust = 2, color = "darkorchid4"
    ) +
    geom_hline(yintercept = inferred_minCounts[tag_y], color = "blue4") +
    annotate("text",
      x = 0, y = inferred_minCounts[tag_y],
      label = paste0(names(inferred_minCounts[tag_y]), ": ", inferred_minCounts[tag_y]),
      vjust = -0.9, hjust = 0.2, color = "blue4"
    ) +
    ggtitle(paste0("Only cells with HTO_maxID == ", tag_x, " or HTO_maxID == ", tag_y))
  scatter_plot3
  file_name3 <- paste(plot_prefix, "scatter_raw_counts", paste0(tag_y, "-highlighted"), tag_x, tag_y, "png", sep = ".")
  print(file_name3)
  ggplot2::ggsave(filename = file_name3, plot = scatter_plot3, width = 19, height = 17, units = "cm")


  # give one plot with all plots
  all_pair_plots <- list(
    scatter_plot_all, scatter_plot_doublet_details, scatter_plot,
    scatter_plot2, scatter_plot3
  )
  my_grid <- cowplot::plot_grid(plotlist = all_pair_plots, ncol = 2)
  file_name_ALL <- paste(opt$output_prefix, "all_paired_scatterplots", tag_x, tag_y, "png", sep = ".")
  print(file_name_ALL)
  ggplot2::ggsave(filename = file_name_ALL, plot = my_grid, width = 50, height = 17 * ceiling(length(all_pair_plots) / 2), units = "cm")
}

##########
# UMAPs
# Calculate UMAP embeddings with a distance matrix
##########
print("Start Calculating UMAP")
colours_hash.IDs <- hue_pal()(length(levels(cellData_htoCounts$HTO_classification)) - 1)
names(colours_hash.IDs) <- levels(cellData_htoCounts$HTO_classification)[levels(cellData_htoCounts$HTO_classification) != "Negative"]
colours_hash.IDs <- c(colours_hash.IDs, "Singlet" = "#b2df8a", "Doublet" = "#1f78b4", "Negative" = "#a6cee3")

cat("UMAP is calculated on raw counts.\n")
umap_input <- hto_counts
umap_input$barcodes <- NULL
umap_input <- umap_input
d <- dist(umap_input)
umap.hto <- umap(d, scale = TRUE, n_neighbors = 3000, pca = 50, a = 2, b = 1.5, fast_sgd = TRUE, ret_nn = T)
umap_coord <- as.data.frame(umap.hto$embedding)
colnames(umap_coord) <- c("umap1", "umap2")
umap_coord$barcodes <- rownames(umap_input)
cellData_htoCounts_extended <- dplyr::inner_join(cellData_htoCounts, umap_coord, by = "barcodes")
# write full cell meta data table into file
filename_all_results <- paste0(opt$output_prefix, ".all_cell_meta_data.tsv")
print(filename_all_results)
write.table(
  x = cellData_htoCounts_extended, file = filename_all_results,
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE
)

# colour setting
cellData_htoCounts_colours <- match(cellData_htoCounts$HTO_maxID, names(colours_hash.IDs))
cluster_maxID <- ggplot(cellData_htoCounts_extended, aes(x = umap1, y = umap2, colour = HTO_maxID)) +
  geom_point(size = 0.005) +
  scale_colour_manual(name = "HTO_maxID", values = colours_hash.IDs[cellData_htoCounts_colours]) +
  theme(
    aspect.ratio = 1,
    legend.key.size = unit(0.8, "line")
  ) +
  # coord_fixed(ratio = 1) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  xlab("UMAP 1") +
  ylab("UMAP 2")
cluster_maxID
cellData_htoCounts_colours <- match(cellData_htoCounts$HTO_classification, names(colours_hash.IDs))
cluster_classification <- ggplot(cellData_htoCounts_extended, aes(x = umap1, y = umap2, color = HTO_classification)) +
  geom_point(size = 0.005) +
  scale_colour_manual(name = "HTO_classification", values = colours_hash.IDs[cellData_htoCounts_colours]) +
  theme(aspect.ratio = 1) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  # coord_fixed(ratio = 1) +
  xlab("UMAP 1") +
  ylab("UMAP 2")
cellData_htoCounts_colours <- match(cellData_htoCounts$HTO_classification.global, names(colours_hash.IDs))
cluster_classification.global <- ggplot(cellData_htoCounts_extended, aes(x = umap1, y = umap2, color = HTO_classification.global)) +
  geom_point(size = 0.005) +
  theme(aspect.ratio = 1) +
  scale_colour_manual(name = "HTO_classification.global", values = colours_hash.IDs[cellData_htoCounts_colours]) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  # coord_fixed(ratio = 1) +
  xlab("UMAP 1") +
  ylab("UMAP 2")
cellData_htoCounts_colours <- match(cellData_htoCounts$hash.ID, names(colours_hash.IDs))
cluster_hashID <- ggplot(cellData_htoCounts_extended, aes(x = umap1, y = umap2, color = hash.ID)) +
  geom_point(size = 0.001) +
  scale_colour_manual(name = "Hash.ID", values = colours_hash.IDs[cellData_htoCounts_colours]) +
  theme(aspect.ratio = 1) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  # coord_fixed(ratio = 1) +
  xlab("UMAP 1") +
  ylab("UMAP 2")

cluster_counts <- ggplot(cellData_htoCounts_extended, aes(x = umap1, y = umap2, colour = log(nCount_HTO, 10))) +
  geom_point(size = rel(0.5)) +
  theme(aspect.ratio = 1) +
  scale_color_gradient(low = "#6a3d9a", high = "#ff7f00") +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  # coord_fixed(ratio = 1) +
  xlab("UMAP 1") +
  ylab("UMAP 2")
ggplot2::ggsave(filename = paste(opt$output_prefix, "cluster_counts.png", sep = "."), plot = cluster_counts, width = 19, height = 17, units = "cm")

cluster_overview <- cluster_classification.global + cluster_classification + cluster_hashID + cluster_maxID
ggplot2::ggsave(filename = paste(opt$output_prefix, "cluster_overview.png", sep = "."), plot = cluster_overview, width = 29, height = 19, units = "cm")


# First, we will remove negative cells from the object. Might not be a good thing to do, thats wy the .subset gets assigned the whole sample currently. Might want to do this later on.
# TODO: Solve this comment?
Idents(my_sample) <- "hash.ID"
# Write out list with negative cells
barcodes_negatives <- names(my_sample$hash.ID[my_sample$hash.ID == "Negative"])
filename_negatives <- paste(opt$output_prefix, "barcodes_negatives.txt", sep = ".")
cat("File with barcodes of Negatives:\n", filename_negatives, "\n")
write.table(barcodes_negatives, filename_negatives, quote = FALSE, row.names = FALSE, col.names = FALSE)
# have Seurat object without Negatives
my_sample.subset <- subset(my_sample, idents = "Negative", invert = TRUE)

# Write out list with Doublets
barcodes_doublets <- names(my_sample$hash.ID[my_sample$hash.ID == "Doublet"])
filename_doublets <- paste(opt$output_prefix, "barcodes_doublets.txt", sep = ".")
cat("File with barcodes of Doublets:\n", filename_doublets, "\n")
write.table(barcodes_doublets, filename_doublets, quote = FALSE, row.names = FALSE, col.names = FALSE)
# have Seurat object without Doublets (or Negatives)
singlets <- my_sample.subset[, my_sample.subset$HTO_classification.global != "Doublet"]

# Write out list with all Singlets
barcodes_singlets <- rownames(singlets[[]])
filename_all_singlets <- paste(opt$output_prefix, "barcodes_singlets.txt", sep = ".")
cat("File with barcodes of Singlets:\n", filename_all_singlets, "\n")
write.table(barcodes_singlets, filename_all_singlets, row.names = FALSE, col.names = FALSE, quote = FALSE)

############
# Write one barcode list for every tag.
###########
cat("\nStart writing out barcode lists.\n\n")
cat("Number all cells: ", length(my_sample$hash.ID), "\n")
cat("Number Negatives: ", length(barcodes_negatives), "\n")
cat("Number Doublets: ", length(barcodes_doublets), "\n")
cat("Number Singlets: ", length(singlets$hash.ID), "\n")
# Check if Singlets, Negatives, and Doublets sum up to the total number of cells.
stopifnot(sum(length(singlets$hash.ID), length(barcodes_doublets), length(barcodes_negatives)) == length(my_sample$hash.ID))

cat("\n\nAll tags in this sample:\n")
print(levels(as.factor(my_sample$HTO_maxID)))
print(str(my_sample$HTO_maxID))

for (tag in levels(as.factor(my_sample$HTO_maxID))) {
  glue("Write into file barcodes of: {tag}")
  # Make sure that we can associate the correct barcode list with a sample.
  sampleName <- sampleTagMap[sampleTagMap$tagName == tag, ]$sample
  print(paste0("Sample name: ", as.character(sampleName)))
  # get subset of all singlets with the current tag
  singlets_current_tag <- singlets[, singlets$hash.ID == tag]
  barcodes_current_tag <- rownames(singlets_current_tag[[]])
  print(paste0("Number of cells with ", tag, ": ", length(barcodes_current_tag)))
  filename_singlets_barcodes <- paste(opt$output_prefix, tag, sampleName, "barcodes_singlets.txt", sep = ".")
  print(filename_singlets_barcodes)
  write.table(barcodes_current_tag, filename_singlets_barcodes, quote = FALSE, row.names = FALSE, col.names = FALSE)
}
