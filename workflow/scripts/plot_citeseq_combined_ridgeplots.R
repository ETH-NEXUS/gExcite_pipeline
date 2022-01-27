## Script to plot combined ridge plot over multiple used tags in one project.
## To do: extend script that it can be used on multiple project having multiple tags.

# Linda Grob, 01.2020

# Example call
# Rscript plotRidgeplots.PBMCmix.20191130.R
# --atypicalRemoved /cluster/project/nexus/phrt/levesque_melanoma_2019/data/20191130/PBMC_mix/analysis_PBMC001/atypical_removed/,/cluster/project/nexus/phrt/levesque_melanoma_2019/data/20191130/PBMC_mix/analysis_PBMC002/atypical_removed/,/cluster/project/nexus/phrt/levesque_melanoma_2019/data/20191130/PBMC_mix/analysis_PBMC003/atypical_removed/
# --analysisADT /cluster/project/nexus/phrt/levesque_melanoma_2019/data/20191130/PBMC_mix/analysis_pooledSample/CellRangerADT/PBMCmix_pooled/outs/filtered_feature_bc_matrix/
# --outfolder /cluster/project/nexus/phrt/levesque_melanoma_2019/data/20191130/PBMC_mix/CompareRidgeplots/ --sampleNames PBMC001,PBMC002,PBMC003


library(cowplot)
library(ggridges)
library(optparse)
library(Seurat)
library(rjson)
library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)
library(limma)
library(viridis)
library(grDevices)
library(patchwork)
library(SingleCellExperiment)

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

option_list <- list(
make_option("--atypicalRemoved", type = "character", help = "Comma separated list with Paths to the atypicalRemoved Folder/RDS Files of the GEX pipline."),
make_option("--analysisADT", type = "character", help = "Path to the pooled ADT analysis folder."),
make_option("--outfolder", type = "character", help = "Path to the outfolder."),
make_option("--sampleNames", type = "character", help = "SampleNames in the order of the atypicalRemoved Input Files.")

)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


outfolder <- opt$outfolder
dir.create(outfolder, showWarnings = FALSE)
outprefix <- paste(outfolder, "CompareRidgeplots", sep = "")
# set general options/parameters
options(stringsAsFactors = FALSE)

# Read in files
#ADT
CellRangerADTfolder <- opt$analysisADT
CellRangerADT <- Read10X(data.dir = CellRangerADTfolder)
# remove "-1" from cell barcodes
colnames(CellRangerADT) <- gsub(pattern = "-1", replacement = "", x = colnames(CellRangerADT))
# shorten AB labels
print("Beautify rownames:")
rownames(x = CellRangerADT) <- gsub(
    pattern = "_TotalSeqC", replacement = "",
    x = rownames(x = CellRangerADT)
)

#GEX: Pre load all rds files, to avoid long runningtime.
my_folders <- as.list(strsplit(opt$atypicalRemoved, ","))[[1]]
my_sampleNames <- as.list(strsplit(opt$sampleNames, ","))[[1]]
my_sceList <- list()
for (tag_folder in 1:length(my_folders)) {
    RDSfile <- list.files(my_folders[tag_folder], pattern = "\\.RDS$")
    my_sceList[[tag_folder]] <- readRDS(paste(my_folders[tag_folder], RDSfile, sep = ""))

}
#Extract Tag spezific ADT subset & Plot subplots
SubmatrixList <- list()

plot_list <- list()
for (i in seq(1, length(rownames(CellRangerADT)))) {
    jointSet <- data.frame(barcode = character(), antibody = character(), value = numeric())
    for (tag_folder in 1:length(my_folders)) {
        my_sce <- my_sceList[[tag_folder]]
        joint.bcs <- intersect(colnames(my_sce), colnames(CellRangerADT))
        CellRangerADTsub <- CellRangerADT[, joint.bcs]
        ABsubset <- as.matrix(CellRangerADTsub[i, ])
        colnames(ABsubset) <- rownames(CellRangerADT)[i]
        ABsubset <- melt(ABsubset)
        colnames(ABsubset) <- c("barcode", "antibody", "value")
        ABsubset$sample <- my_sampleNames[[tag_folder]]
        #ABsubset$value[ABsubset$value==0] <- NA
        # add plus one (+1) to the counts before taking log()
        ABsubset$value <- ABsubset$value + 1
        ABsubset$ln <- log(ABsubset$value)
        jointSet <- rbind(jointSet, ABsubset)
    }
    s <- ggplot(jointSet, aes(x = ln, y = sample, fill = sample, height = ..ndensity..)) +
          geom_density_ridges(alpha = 0.3) +
          scale_x_continuous(limits = c(0, max(jointSet$ln) + 1)) +
          #geom_vline(xintercept = log(threshold)) +
          theme_ridges() +
          ggtitle(rownames(CellRangerADT)[i]) +
          theme(legend.position = "none",
                axis.title.y = element_blank(),
                plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"),
                panel.background = element_rect(fill = "white"),
                plot.background = element_rect(fill = "white"))
    s
    plot_list[[i]] <- s
}

# Save plots in combined plots on 4x4 grids.
numberAntibodies <- length(rownames(CellRangerADT))
print(seq(1, 16 * ceiling(numberAntibodies / 16), 16))
print(seq(1, numberAntibodies, 16))
for (i in seq(1, numberAntibodies, 16)) {
        print(i)
        plot_sub <- plot_list[i:min(i + 15, numberAntibodies)]
        plot_group <-  plot_grid(plotlist = plot_sub, ncol = 4)
        if (i + 15 > numberAntibodies) {
        ggplot2::ggsave(filename = paste(outprefix, i, "to", numberAntibodies, "ridgeplot.png", sep = "_"),
                        width = 21, height = ((3.5 * ceiling((length(plot_sub)) / 4))), plot = plot_group)
        } else {
         ggplot2::ggsave(filename = paste(outprefix, i, "to", i + 15, "ridgeplot.png", sep = "_"),
                         width = 21, height = 14, plot = plot_group)
}
}
