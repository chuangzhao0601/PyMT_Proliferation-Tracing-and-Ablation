# ------------------------------------------------------------------
# Integration and Dimensionality Reduction of Non-Immune and Epithelial Cells in the PyMT Mouse Model
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
##Nonimmune cells and Epithelial

#create Seurat objects
samples <- c("A0305","A0306","A0307","A0308")
seurat_list <- lapply(samples, function(sample){
  cur_data <- Read10X(data.dir = paste0("../PyMT/download_data/",sample))
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
nonimmune <- merged_seurat

# clean up
rm(seurat_list)
rm(merged_seurat)

#remove the contamination of ambient RNA using DecontX
nonimmune_sce <- SingleCellExperiment(list(counts = nonimmune@assays$RNA@counts))
nonimmune_sce <- decontX(nonimmune_sce)
nonimmune[["decontXcounts"]] <- CreateAssayObject(counts = nonimmune_sce@assays@data$decontXcounts)
DefaultAssay(nonimmune) <- "decontXcounts"
nonimmune_sce <- SingleCellExperiment(list(counts = nonimmune@assays$decontXcounts@counts))
nonimmune_sce <- decontX(nonimmune_sce)
nonimmune[["decontXcounts"]] <- CreateAssayObject(counts = nonimmune_sce@assays@data$decontXcounts)
DefaultAssay(nonimmune) <- "decontXcounts"


#label the intervene
nonimmune@meta.data$intervene <- NA
nonimmune@meta.data$intervene[nonimmune@meta.data$orig.ident == "A0305"] <- "TAM"
nonimmune@meta.data$intervene[nonimmune@meta.data$orig.ident == "A0306"] <- "TAM"
nonimmune@meta.data$intervene[nonimmune@meta.data$orig.ident == "A0307"] <- "TAM_DT"
nonimmune@meta.data$intervene[nonimmune@meta.data$orig.ident == "A0308"] <- "TAM_DT"

#Quality control
nonimmune[["percent.mt"]] <- PercentageFeatureSet(nonimmune, pattern = "^MT-|^mt-")
nonimmune[["percent.rb"]] <- PercentageFeatureSet(nonimmune, pattern = "^RP[SL]|^Rp[sl]")
VlnPlot(object = nonimmune, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"),group.by="orig.ident", ncol = 3)
nonimmune <- subset(nonimmune, subset =nFeature_RNA > 200 & nFeature_RNA <= quantile(nFeature_RNA, 0.98) 
                    & nCount_RNA <= quantile(nCount_RNA, 0.98)
                    & percent.mt < 10)

# Normalize the data
nonimmune <- NormalizeData(nonimmune, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable genes for PCA 
nonimmune <- FindVariableFeatures(object = nonimmune,selection.method = "vst", nfeatures = 2000)
# Scale the data
all.gene <- rownames(nonimmune)
nonimmune <- ScaleData(nonimmune,features = all.gene)

# Perform linear dimensional reduction
nonimmune <- RunPCA(object=nonimmune,features = VariableFeatures(object = nonimmune))
DimPlot(nonimmune, reduction = "pca")
ElbowPlot(nonimmune,n=30)

# Cluster the cells
nonimmune <- FindNeighbors(nonimmune, reduction = "pca", dims = 1:30)
nonimmune <- FindClusters(nonimmune,resolution = 1.5, algorithm = 1)

# Run non-linear dimensional reduction (tSNE and uMAP)
nonimmune <- RunTSNE(object=nonimmune,dims.use=1:30,do.fast=TRUE,check_duplicates = FALSE)
nonimmune <- RunUMAP(nonimmune, reduction = "pca", dims = 1:30)

DimPlot(nonimmune, reduction = "umap",label = TRUE)

# Find markers 
markers <- FindAllMarkers(object = nonimmune, 
                          only.pos = TRUE,
                          logfc.threshold= 0.25,
                          min.pct = 0.25)

#write the markers of each cluster
write.table(markers,file="nonimmune_marker.tsv",sep="\t")

Idents(nonimmune) <- "seurat_clusters"
new.cluster.ids <- c("0"="Epithelial", 
                     "1"="Uncharacterized", 
                     "2"="Epithelial", 
                     "3"="Epithelial", 
                     "4"="Uncharacterized", 
                     "5"="Epithelial", 
                     "6"="Immune", 
                     "7"="Epithelial", 
                     "8"="Immune", 
                     "9"="Fibroblast", 
                     "10"="Epithelial", 
                     "11"="Uncharacterized", 
                     "12"="Uncharacterized", 
                     "13"="Epithelial", 
                     "14"="Epithelial", 
                     "15"="Epithelial", 
                     "16"="Epithelial", 
                     "17"="Uncharacterized",
                     "18"="Immune", 
                     "19"="Immune", 
                     "20"="Endothelial",
                     "21"="Pericyte", 
                     "22"="Myoepithelial", 
                     "23"="Myoepithelial",
                     "24"="Uncharacterized",
                     "25"="Fibroblast",
                     "26"="Epithelial",
                     "27"="Immune"
)
# Add cell type to the metadata as a new column
nonimmune <- RenameIdents(nonimmune, new.cluster.ids)
nonimmune$celltype <- nonimmune@active.ident

#Extract the real nonimmune cells
nonimmune <- subset(nonimmune,celltype%in%c("Epithelial","Fibroblast","Myoepithelial","Endothelial","Pericyte"))

#Save the RDS
saveRDS(nonimmune,"nonimmune_annotion.rds")

#Annotion the Epithelial####
nonimmune <- readRDS("../nonimmune_annotion.rds")
unique(nonimmune$celltype)
Epithelial <- subset(nonimmune,celltype%in%c("Epithelial"))
rm(nonimmune)
# Find variable genes for PCA 
Epithelial <- FindVariableFeatures(object = Epithelial,selection.method = "vst", nfeatures = 2000)

# Scale the data
all.gene <- rownames(Epithelial)
Epithelial <- ScaleData(Epithelial,features = all.gene)

# Perform linear dimensional reduction
Epithelial <- RunPCA(object=Epithelial,features = VariableFeatures(object = Epithelial))
DimPlot(Epithelial, reduction = "pca")
ElbowPlot(Epithelial,n=30)

# Cluster the cells
Epithelial <- FindNeighbors(Epithelial, reduction = "pca", dims = 1:30)
Epithelial <- FindClusters(Epithelial,resolution = 1, algorithm = 1)

# Run non-linear dimensional reduction (tSNE and uMAP)
Epithelial <- RunTSNE(object=Epithelial,dims.use=1:30,do.fast=TRUE,check_duplicates = FALSE)
Epithelial <- RunUMAP(Epithelial, reduction = "pca", dims = 1:30,min.dist = 1,n.neighbors = 40L)
DimPlot(Epithelial, reduction = "umap",label = TRUE)

# Find markers 
markers <- FindAllMarkers(object = Epithelial, 
                          only.pos = TRUE,
                          logfc.threshold= 0.25,
                          min.pct = 0.25)

#write the markers of each cluster
write.table(markers,file="Epithelial_dims30_resolution1.tsv",sep="\t")

#annotion
Idents(Epithelial) <- "seurat_clusters"
new.cluster.ids <- c("0"="Cancer_cells_c1_Spp1", 
                     "1"="Cancer_cells_c2_Tspan8", 
                     "2"="Cancer_cells_c3_Hspa1a", 
                     "3"="Cancer_cells_c4_Cxcl1",
                     "4"="Cancer_cells_c5_Lalba", 
                     "5"="Cancer_cells_c6_Pde4d", 
                     "6"="Cancer_cells_c7_Myl9", 
                     "7"="Cancer_cells_c8_Ntrk2",
                     "8"="Cancer_cells_c9_Top2a", 
                     "9"="Cancer_cells_c10_Ly6a", 
                     "10"="Cancer_cells_c11_Ifit3",
                     "11"="Cancer_cells_c12_Aldh1a3", 
                     "12"="Cancer_cells_c13_Prlr", 
                     "13"="Cancer_cells_c14_Csn1s1"
)
# Add cell type to the metadata as a new column
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)
Epithelial$subcelltype1 <- Epithelial@active.ident

#Save the RDS
saveRDS(Epithelial,"PyMT_Epithelial.rds")
