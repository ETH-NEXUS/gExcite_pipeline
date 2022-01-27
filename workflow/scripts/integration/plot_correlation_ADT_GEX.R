# Script to plot correlation between ADT and GEX expression
# Produces a table and a plot
# Input is an integrated seurat object and a threshold as well as a lookup table.


# ToDo 
# Remove unnecessary libs
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
  library(reshape2)
  library(dplyr)
  library(tidyverse)
})

# parse command line arguments
 option_list <- list(
   make_option("--cohort_object", type = "character", help = "Integrated Seurat Object"),
   make_option("--lookup", type = "character", help = "Full path to a lookup table matching protein to antibody names."),
   make_option("--thresholds", type = "character", help = "Full path to a table listing adt thresholds for every sample included in the annalysis."),
   make_option("--outdir", type = "character", help = "Full path to output directory."),
   make_option("--outName", type = "character", help = "Prefix name of output files.")
 )
 opt_parser <- OptionParser(option_list = option_list)
 opt <- parse_args(opt_parser)


# Read cohort object
seurat_integrated <- readRDS(opt$cohort_object)
print(seurat_integrated)

# Read in lookup table

lookup <-read.csv(opt$lookup, sep = "\t", stringsAsFactors = FALSE)
lookup$Antibody <- sub("_", "-", lookup$Antibody)
print(lookup)

# Read in thresholds
thresholds <- read.csv(opt$thresholds, sep= "\t", stringsAsFactors = FALSE)
thresholds$Antibody<- sub("_", "-", thresholds$Antibody)
rownames(thresholds) <- toupper(thresholds$Antibody)
print(tail(thresholds,20))

# Make new assay with threshold
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

my_new_assay <- my_merge[,1:107]
rownames(my_new_assay) <- my_new_assay$Row.names
my_new_assay$Row.names <- NULL
seurat_integrated@assays$adt_thresholds <-  CreateAssayObject(counts = t(my_new_assay))
correlation <- data.frame(matrix(ncol = 6, nrow = 0))

#provide column names
colnames(correlation) <- c('ab', 'gene', 'both_expressed', 'both_not_expressed', 'only_ab', 'only_gex')

for (ab in seq(1,nrow(lookup))){
  print(ab)
  if (lookup$Antibody[ab] %in% c("hashtag1","hashtag2","hashtag3","hashtag4")){
    print(lookup$Antibody[ab])
    print("hashtag")
    next
  }
  gene <- lookup$Gene[ab]
  print(gene)
  gene_expressed <- TRUE
  if (!gene %in% rownames(seurat_integrated@assays$SCT@counts)){
    gene_expressed <- FALSE
  }
  both_expressed <- 0
  both_not_expressed <- 0
  only_ab <- 0
  only_gex <- 0
  numberCells <- ncol(seurat_integrated@assays$adt_thresholds@counts)
  
  pos_for_ab <- sum(seurat_integrated@assays$adt_thresholds@counts[paste0("ADT-",lookup$Antibody[ab]),]>0)
  pos_barcodes_ab <-  seurat_integrated@assays$adt_thresholds@counts[paste0("ADT-",lookup$Antibody[ab]),]>0
  pos_barcodes_ab <- names(pos_barcodes_ab[pos_barcodes_ab==TRUE])
  if (gene_expressed){
  pos_for_gex <- sum(seurat_integrated@assays$SCT@counts[gene,]>0)
  pos_barcodes_gex <- seurat_integrated@assays$SCT@counts[gene,]>0
  pos_barcodes_gex <- names(pos_barcodes_gex[pos_barcodes_gex==TRUE])
  both_expressed <- length(intersect(pos_barcodes_gex, pos_barcodes_ab))
  only_ab <- length(setdiff(pos_barcodes_ab, pos_barcodes_gex))
  only_gex <- length(setdiff(pos_barcodes_gex, pos_barcodes_ab))
  both_not_expressed <- length(setdiff(rownames(seurat_integrated@meta.data), union(pos_barcodes_ab, pos_barcodes_gex)))
  }
  else{
    both_expressed <- 0
    only_gex <- 0
    only_ab <- length(pos_barcodes_ab)
    both_not_expressed <- length(setdiff(rownames(seurat_integrated@meta.data), pos_barcodes_ab))
    
  }
  print(paste0("Expressed in both: ", both_expressed))
  print(paste0("Expressed in none: ", both_not_expressed))
  print(paste0("Expressed in only ab: ", only_ab))
  print(paste0("Expressed in only gex: ", only_gex))
  print(paste0("Total: ", both_expressed+both_not_expressed+ only_ab+only_gex))
  correlation<-rbind(correlation,  c(lookup$Antibody[ab], gene, both_expressed, both_not_expressed, only_ab, only_gex))
}

#provide column names
colnames(correlation) <- c('ab', 'gene', 'c_both_expressed', 'd_both_not_expressed', 'b_only_ab', 'a_only_gex')

correlation_long <- correlation
correlation_long$identifier <- paste0(correlation_long$ab,"/", correlation_long$gene)
for (item in seq(1,nrow(correlation_long))){
  if (toupper(correlation_long$ab[item]) == toupper(correlation_long$gene[item])){
    correlation_long$identifier[item] <- toupper(correlation_long$ab[item])
  } 
}


correlation_long$ab <- NULL
correlation_long$gene <- NULL
correlation_long$d_both_not_expressed <- as.integer(correlation$d_both_not_expressed)


correlation_long <- melt(correlation_long, id.vars = "identifier")
correlation_long$value <- as.integer(correlation_long$value)
correlation_long$identifier <- as.factor(correlation_long$identifier)
correlation_ordered <- arrange(correlation_long[correlation_long$variable=="d_both_not_expressed",], value)
correlation_long$variable <- relevel(correlation_long$variable,ref= c("d_both_not_expressed"))
# Reorder following the value of another column:
correlation_plot <- ggplot(    correlation_long[order(correlation_long$variable, decreasing = TRUE),],            # Stacked barplot using ggplot2
       aes(x = identifier,
           y = value,
           fill = variable)) +
  geom_bar(position="fill",stat = "identity") +
  scale_x_discrete(limits = correlation_ordered$identifier) +
  labs(x ="Antibody / Gene", y = "Percentage", fill = "Correlation") +
  scale_fill_manual(labels = c("Gene and AB not expressed", "Gene and AB expressed", "Only AB expressed", "Only Gene expressed"),values = c("#1B9E77", "#D95F02" ,"#7570B3" , "#66A61E"))+
theme(axis.text.x = element_text(angle = 90,size = 15, vjust = 0.5, hjust=1), panel.background = element_blank(),
      axis.title.x = element_text(size = 25),
      axis.title.y = element_text(size = 25),
      axis.text.y= element_text(size = 25),
      legend.text = element_text(size = 25),
      legend.title = element_text(size = 25))
ggsave(filename = paste0(opt$outdir, opt$outName, "_correlation.png"),
       width = nrow(lookup)*0.7, height = 20, dpi = 600, units = "cm")


colnames(correlation) <- c("Antibody","Gene","Gene and AB expressed", "Gene and AB not expressed", "Only AB expressed", "Only Gene expressed")
write.table(correlation, sep="\t", quote=FALSE,file = paste0(opt$outdir, opt$outName, "_correlation_GEX_AB.txt"), row.names = FALSE)
