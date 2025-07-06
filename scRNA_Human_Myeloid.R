# ------------------------------------------------------------------
# Reanalysis of Myeloid Cells from Human Breast Cancer scRNA-seq Data
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

#Read the RDS
bc2021 <- readRDS("/home/shixi7/zhaochuang/project/pymt/PyMT/human/BC_2021_all_anno.rds")
DimPlot(bc2021, reduction = "umap",label = TRUE,group.by = "cellType")
unique(bc2021$cellType)

#Select the Myeloid Cell
Myeloid_cell <- subset(bc2021,cellType%in%c("Myeloid_cell"))
rm(bc2021)

#Quality control
Myeloid_cell[["percent.mt"]] <- PercentageFeatureSet(Myeloid_cell, pattern = "^MT-|^mt-")
Myeloid_cell[["percent.rb"]] <- PercentageFeatureSet(Myeloid_cell, pattern = "^RP[SL]|^Rp[sl]")
VlnPlot(object = Myeloid_cell, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"),group.by="orig.ident", ncol = 3)
Myeloid_cell <- subset(Myeloid_cell, subset =nFeature_RNA > 200 & nFeature_RNA <6000
                       &nCount_RNA >400)

# Find variable genes for PCA 
Myeloid_cell <- FindVariableFeatures(object = Myeloid_cell,selection.method = "vst", nfeatures = 2000)
# Scale the data
all.gene <- rownames(Myeloid_cell)
Myeloid_cell <- ScaleData(Myeloid_cell,features = all.gene)

# Perform linear dimensional reduction
Myeloid_cell <- RunPCA(object=Myeloid_cell,features = VariableFeatures(object = Myeloid_cell))
DimPlot(Myeloid_cell, reduction = "pca")
ElbowPlot(Myeloid_cell,n=30)

# Cluster the cells
Myeloid_cell <- FindNeighbors(Myeloid_cell, reduction = "pca", dims = 1:30)
Myeloid_cell <- FindClusters(Myeloid_cell,resolution = 1, algorithm = 1)

# Run non-linear dimensional reduction (tSNE and uMAP)
Myeloid_cell <- RunTSNE(object=Myeloid_cell,dims.use=1:30,do.fast=TRUE,check_duplicates = FALSE)
Myeloid_cell <- RunUMAP(Myeloid_cell, reduction = "pca", dims = 1:30)
DimPlot(Myeloid_cell, reduction = "umap",label = TRUE)

# Find markers 
markers <- FindAllMarkers(object = Myeloid_cell, 
                          only.pos = TRUE,
                          logfc.threshold= 0.25,
                          min.pct = 0.25)

#write the markers of each cluster
write.table(markers,file="BC2021_marker.tsv",sep="\t")

#Annotion
Idents(Myeloid_cell) <- "seurat_clusters"
new.cluster.ids <- c("0"="MPS", 
                     "1"="Uncharacterized", 
                     "2"="MPS", 
                     "3"="MPS", 
                     "4"="Uncharacterized", 
                     "5"="MPS", 
                     "6"="MPS", 
                     "7"="MPS", 
                     "8"="Uncharacterized", 
                     "9"="MPS", 
                     "10"="Uncharacterized", 
                     "11"="Uncharacterized", 
                     "12"="MPS", 
                     "13"="MPS", 
                     "14"="T", 
                     "15"="Uncharacterized", 
                     "16"="MPS", 
                     "17"="Uncharacterized",
                     "18"="Uncharacterized", 
                     "19"="Uncharacterized", 
                     "20"="Uncharacterized",
                     "21"="Uncharacterized", 
                     "22"="Uncharacterized", 
                     "23"="Uncharacterized",
                     "24"="Uncharacterized",
                     "25"="Uncharacterized"
)
# Add cell type to the metadata as a new column
Myeloid_cell <- RenameIdents(Myeloid_cell, new.cluster.ids)
Myeloid_cell$celltype <- Myeloid_cell@active.ident

#Select the MPS
MPS <- subset(Myeloid_cell,celltype%in%c("MPS"))

# Find variable genes for PCA 
MPS <- FindVariableFeatures(object = MPS,selection.method = "vst", nfeatures = 2000)
# Scale the data
all.gene <- rownames(MPS)
MPS <- ScaleData(MPS,features = all.gene)

# Perform linear dimensional reduction
MPS <- RunPCA(object=MPS,features = VariableFeatures(object = MPS))
DimPlot(MPS, reduction = "pca")
ElbowPlot(MPS,n=30)

# Cluster the cells
MPS <- FindNeighbors(MPS, reduction = "pca", dims = 1:30)
MPS <- FindClusters(MPS,resolution = 1.5, algorithm = 1)

# Run non-linear dimensional reduction (tSNE and uMAP)
MPS <- RunTSNE(object=MPS,dims.use=1:30,do.fast=TRUE,check_duplicates = FALSE)
MPS <- RunUMAP(MPS, reduction = "pca", dims = 1:30)
DimPlot(MPS, reduction = "umap",label = TRUE)

# Find markers 
markers <- FindAllMarkers(object = MPS, 
                          only.pos = TRUE,
                          logfc.threshold= 0.25,
                          min.pct = 0.25)

#write the markers of each cluster
write.table(markers,file="MPS_marker.tsv",sep="\t")

#Annotion
Idents(MPS) <- "seurat_clusters"
new.cluster.ids <- c("0"="Mac_Other", 
                     "1"="Mac_Other", 
                     "2"="Mac_Other", 
                     "3"="Mac_Other", 
                     "4"="Mac_Other", 
                     "5"="Mono", 
                     "6"="Mono", 
                     "7"="Mac_SPP1_VEGFA", 
                     "8"="Mac_Other", 
                     "9"="Mac_Other", 
                     "10"="Mac_Other", 
                     "11"="Mac_Other", 
                     "12"="Mac_Other", 
                     "13"="Mac_SPP1_VEGFA", 
                     "14"="Mac_Other", 
                     "15"="Mac_Other", 
                     "16"="Mac_SPP1_VEGFA", 
                     "17"="Mac_Other",
                     "18"="Mono", 
                     "19"="Mac_SPP1_VEGFA", 
                     "20"="Mac_Other",
                     "21"="Mac_Other", 
                     "22"="Mac_Other"
)
# Add cell type to the metadata as a new column
MPS <- RenameIdents(MPS, new.cluster.ids)
MPS$subcelltype1 <- MPS@active.ident

#Select the required cells
MPS <- subset(MPS,expansion%in%c("E","NE"))
MPS <- subset(MPS,BC_type%in%c("ER+","TNBC"))
MPS <- subset(MPS,subcelltype1%in%c("Mac_SPP1_VEGFA","Mac_Other"))

#Create a new column
MPS$BC_type_expansion <- paste0(MPS$BC_type,sep="_",MPS$expansion)

#Save the RDS
saveRDS(MPS,"Human_bc2021_MPS.rds")