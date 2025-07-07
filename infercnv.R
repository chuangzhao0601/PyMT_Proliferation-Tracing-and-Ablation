# ------------------------------------------------------------------
# InferCNV analysis on tumor and immune cells (mouse data)
# ------------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(infercnv)
library(AnnoProbe)
library(magrittr)
library(scales)
library(ggpubr)
library(tidyr)

# Load Seurat object
all <- readRDS("../PyMT_all.rds")

# Extract reference immune cells (Macrophages)
Immune_cells <- rownames(all@meta.data)[all$celltype == 'Macrophage']
Immune_cells_Mat <- as.data.frame(GetAssayData(subset(all, cells = Immune_cells), slot = 'counts'))

# Extract Ly6a+ tumor subcluster
ly6a <- rownames(all@meta.data)[all$subcelltype1 == 'Cancer_cells_c10_Ly6a']
ly6a_Mat <- as.data.frame(GetAssayData(subset(all, cells = ly6a), slot = 'counts'))

# Extract other tumor cells (excluding Ly6a)
other <- rownames(all@meta.data)[all$subcelltype1 != 'Cancer_cells_c10_Ly6a']
other_Mat <- as.data.frame(GetAssayData(subset(all, cells = other), slot = 'counts'))

# Combine all matrices
dat <- cbind(Immune_cells_Mat, ly6a_Mat, other_Mat)

# Create annotation file (cell barcode to group label)
groupinfo <- data.frame(
  v1 = colnames(dat),
  v2 = c(
    rep('Immune_cells', ncol(Immune_cells_Mat)),
    rep('Cancer_cells_ly6a', ncol(ly6a_Mat)),
    rep('Cancer_cells_other', ncol(other_Mat))
  )
)

# Get gene positions (SYMBOL to genomic coordinates)
geneInfor <- annoGene(rownames(dat), "SYMBOL", "mouse")
geneInfor <- geneInfor[order(geneInfor$chr, geneInfor$start), c(1, 4:6)]
geneInfor <- geneInfor[!duplicated(geneInfor[,1]), ]

# Filter gene matrix by those with known coordinates
dat <- dat[rownames(dat) %in% geneInfor[,1], ]

# Set output dir
setwd("/home/shixi7/zhaochuang/project/pymt/PyMT/normal-annotion/infercnv2/")

# Save required input files
write.table(dat, file = "expFile.txt", sep = "\t", quote = FALSE)
write.table(groupinfo, file = "groupFiles.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(geneInfor, file = "geneFile.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# ------------------------------------------------------------------
# Create InferCNV object
# ------------------------------------------------------------------

expFile <- "expFile.txt"
groupFiles <- "groupFiles.txt"
geneFile <- "geneFile.txt"

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = expFile,
  annotations_file = groupFiles,
  gene_order_file = geneFile,
  ref_group_names = "Immune_cells"
)
saveRDS(infercnv_obj, "infercnv.RDS")

# ------------------------------------------------------------------
# Run InferCNV
# ------------------------------------------------------------------

pro <- readRDS("infercnv.RDS")
infercnv_obj3 <- infercnv::run(
  pro,
  cutoff = 0.1,                         # min average read counts per gene
  denoise = TRUE,                       # remove noise
  sd_amplifier = 1.5,                   # amplification threshold
  HMM = FALSE,                          # skip HMM inference
  k_obs_groups = 1,
  cluster_by_groups = TRUE,
  analysis_mode = "subclusters",
  tumor_subcluster_partition_method = "qnorm",
  scale_data = TRUE,
  write_expr_matrix = TRUE,
  num_threads = 10,
  out_dir = "CNVout_ceshi"
)
saveRDS(infercnv_obj3, "ceshi_infercnv.RDS")
