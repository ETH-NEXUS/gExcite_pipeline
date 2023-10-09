##################################################
#### File name: dsb_normalize_adt.R
#### Author: Anne Bertolini
#### October 2023
#### R Version: 4
##################################################

library(dsb)
library(Seurat)
library(ggplot2)

# command line arguments are parsed
option_list <- list(
  make_option("--adt_raw_dir", type = "character", help = "Path to directory from cellranger, raw output, \"raw_feature_bc_matrix\"."),
  make_option("--adt_filtered_dir", type = "character", help = "Path to directory from cellranger, output with called cells, \"filtered_feature_bc_matrix\""),
  make_option("--outdir", type = "character", help = "Path to the output directory."),
  make_option("--sample", type = "character", help = "Sample identifier.")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# give out date, time and session info
cat("\n\n")
print(Sys.time())
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n")
cat("\nInput files:\n\n")
print(opt)
cat("\n\n")


# read raw ADT data using the Seurat function "Read10X"
raw_adt = Seurat::Read10X(opt$adt_raw_dir, strip.suffix = TRUE)
cells_adt = Seurat::Read10X(opt$adt_filtered_dir, strip.suffix = TRUE)

# define cell-containing barcodes and separate cells and empty drops
stained_cells_adt = colnames(cells_adt)
background_adt = setdiff(colnames(raw_adt), stained_cells_adt)

# create metadata of droplet QC stats
md = data.frame(
  prot.size = log10(Matrix::colSums(raw_adt))
)

# add indicator for barcodes Cell Ranger called as cells
md$drop.class = ifelse(rownames(md) %in% stained_cells_adt, 'cell', 'background')

# remove barcodes with no evidence of capture in the experiment
md = md[md$prot.size > 0, ]
# remove barcodes from matrix
raw_adt <- raw_adt[, colnames(raw_adt) %in% rownames(md)]

# plot general QC
plot_density <- ggplot2::ggplot(md, aes(x = seq_along(prot.size), y = prot.size )) +
  theme_bw() +
  geom_bin2d(bins = 300) +
  scale_fill_viridis_c(option = "C") +
  facet_wrap(~drop.class)
filename <- paste0(opt$outdir, opt$sample, ".plot_adt_counts.png")
ggsave(filename = filename, plot = plot_density)

# have matrix of background and cells
cells.adt.mtx = as.matrix(raw_adt[, colnames(raw_adt) %in% stained_cells_adt])
background.adt.mtx = as.matrix(raw_adt[, colnames(raw_adt) %in% background_adt])


# give out antibodies with lowest counts in cells, and highest
cat("\n\n\nAntibodies with the lowest maximum counts in sample:\n\n")
pm = sort(apply(cells.adt.mtx, 1, max))
print(head(pm))
cat("\n\n\nAntibodies with the highest maximum counts in sample:\n\n")
print(tail(pm))



# normalize with dsb (denoising is only possible with isoform controls) with
cells.dsb.norm = dsb::DSBNormalizeProtein(
  cell_protein_matrix = cells.adt.mtx,
  empty_drop_matrix = background.adt.mtx,
  denoise.counts = FALSE,
  return.stats = TRUE
)

# write dsb-normalized ADT counts
filename <- paste0(opt$outdir, opt$sample, ".dsb_normalize_adt.RDS")
saveRDS(cells.dsb.norm, file = filename)
