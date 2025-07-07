# ------------------------------------------------------------------
# Integration and Dimensionality Reduction of Immune Cells in the PyMT Mouse Model
# ------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(reshape2)
library(dplyr)
library(grid)
library(tidyverse)
library(biomaRt)
library(celda)
library(SingleCellExperiment)
library(scater)
##Immune cells

#create Seurat objects
samples <- c("A0305","A0306","A0307","A0308")
seurat_list <- lapply(samples, function(sample){
  cur_data <- Read10X(data.dir = paste0("../download_CD45/",sample))
  cur_seurat <- CreateSeuratObject(
    counts = cur_data,
    min.cells=3,
    min.features=200,
    project=sample
  )
  cur_seurat$SampleID <- sample
  return(cur_seurat)
})

# merge Seurat objects
merged_seurat <- merge(x=seurat_list[[1]], y=seurat_list[2:length(seurat_list)])
immune <- merged_seurat

# clean up
rm(seurat_list)
rm(merged_seurat)

#remove the contamination of ambient RNA using DecontX
immune_sce <- SingleCellExperiment(list(counts = immune@assays$RNA@counts))
immune_sce <- decontX(immune_sce)
immune[["decontXcounts"]] <- CreateAssayObject(counts = immune_sce@assays@data$decontXcounts)
DefaultAssay(immune) <- "decontXcounts"
immune_sce <- SingleCellExperiment(list(counts = immune@assays$decontXcounts@counts))
immune_sce <- decontX(immune_sce)
immune[["decontXcounts"]] <- CreateAssayObject(counts = immune_sce@assays@data$decontXcounts)
DefaultAssay(immune) <- "decontXcounts"

#label the intervene
immune@meta.data$intervene <- NA
immune@meta.data$intervene[immune@meta.data$orig.ident == "A0305"] <- "TAM"
immune@meta.data$intervene[immune@meta.data$orig.ident == "A0306"] <- "TAM"
immune@meta.data$intervene[immune@meta.data$orig.ident == "A0307"] <- "TAM_DT"
immune@meta.data$intervene[immune@meta.data$orig.ident == "A0308"] <- "TAM_DT"

#Quality control
immune[["percent.mt"]] <- PercentageFeatureSet(immune, pattern = "^MT-|^mt-")
immune[["percent.rb"]] <- PercentageFeatureSet(immune, pattern = "^RP[SL]|^Rp[sl]")
VlnPlot(object = immune, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"),group.by="orig.ident", ncol = 3)
immune <- subset(immune, subset =nFeature_RNA > 200 & nFeature_RNA <= quantile(nFeature_RNA, 0.98) 
                 & nCount_RNA <= quantile(nCount_RNA, 0.98) 
                 & percent.mt < 10)

# Normalize the data
immune <- NormalizeData(immune, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable genes for PCA 
immune <- FindVariableFeatures(object = immune,selection.method = "vst", nfeatures = 2000)

# Scale the data
all.gene <- rownames(immune)
immune <- ScaleData(immune,features = all.gene)

# Perform linear dimensional reduction
immune <- RunPCA(object=immune,features = VariableFeatures(object = immune))
DimPlot(immune, reduction = "pca")
ElbowPlot(immune,n=30)

# Cluster the cells
immune <- FindNeighbors(immune, reduction = "pca", dims = 1:30)
immune <- FindClusters(immune,resolution = 1.5, algorithm = 1)

# Run non-linear dimensional reduction (tSNE and uMAP)
immune <- RunTSNE(object=immune,dims.use=1:30,do.fast=TRUE,check_duplicates = FALSE)
immune <- RunUMAP(immune, reduction = "pca", dims = 1:30)

DimPlot(immune, reduction = "umap",label = TRUE)

Idents(immune) <- "seurat_clusters"
# Find markers 
markers <- FindAllMarkers(object = immune, 
                          only.pos = TRUE,
                          logfc.threshold= 0.25,
                          min.pct = 0.25)

#write the markers of each cluster
write.table(markers,file="immune_dims30_resolution1.5.tsv",sep="\t")

new.cluster.ids <- c("0"="Tcell", 
                     "1"="Tcell", 
                     "2"="Tcell", 
                     "3"="NK", 
                     "4"="Macrophage", 
                     "5"="Tcell", 
                     "6"="Tcell", 
                     "7"="Tcell", 
                     "8"="Tcell", 
                     "9"="DC", 
                     "10"="Tcell", 
                     "11"="Nonimmune", 
                     "12"="Neutrophil", 
                     "13"="Nonimmune", 
                     "14"="Tcell", 
                     "15"="Neutrophil", 
                     "16"="Monocyte", 
                     "17"="Macrophage",
                     "18"="Bcell", 
                     "19"="Macrophage", 
                     "20"="DC",
                     "21"="Tcell", 
                     "22"="DC", 
                     "23"="Nonimmune",
                     "24"="Nonimmune"
)
# Add cell type to the metadata as a new column
immune <- RenameIdents(immune, new.cluster.ids)
immune$celltype <- immune@active.ident

#Extract the real immune cells
immune <- subset(immune,celltype!="Nonimmune")

#Save the RDS
saveRDS(immune,"immune_annotion.rds")
