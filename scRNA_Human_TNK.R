# ------------------------------------------------------------------
# Reanalysis of TNK Cells from Human Breast Cancer scRNA-seq Data
# ------------------------------------------------------------------
library(Seurat)       
library(harmony)      
library(ggplot2)       
library(dplyr)        
library(patchwork)
library(readr)
library(RColorBrewer)
library(future)
library(cowplot)

all <- readRDS("/hwdata/home/zcshixi7/project/pymt/analyze/Tcell/harmony/Tcell.rds")

all$batch <- paste(all$timepoint,all$patient_id,sep = "_")

# Normalize the data
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable genes for PCA 
all <- FindVariableFeatures(object = all,selection.method = "vst", nfeatures = 2000)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

all <- CellCycleScoring(all, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
all <- ScaleData(all, vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt","S.Score","G2M.Score"), features = rownames(all))

# Perform linear dimensional reduction
all <- RunPCA(object=all,features = VariableFeatures(object = all))
DimPlot(all, reduction = "pca")
ElbowPlot(all,n=50)

# Selecting principal components that explain 90% of the total variance
ElbowPlot(all, reduction = "pca", ndims = 50)
xx <- cumsum(all[["pca"]]@stdev^2)
xx <- xx / max(xx)
which(xx > 0.90)

# Integrated using Harmony
library(harmony)
all <- RunHarmony(all,group.by="batch",plot_convergence=TRUE)

# Cluster the cells
all <- FindNeighbors(all, reduction = "harmony", dims = 1:35)
resolutions <- seq(1.0, 2.0, by = 0.2)
all <- FindClusters(all,resolution = resolutions, algorithm = 1, verbose = FALSE)

# Run non-linear dimensional reduction (tSNE and uMAP)
all <- RunTSNE(object=all,dims.use=1:35,do.fast=TRUE,check_duplicates = FALSE,reduction = "harmony")
all <- RunUMAP(all, reduction = "harmony", dims = 1:35)

saveRDS(all,"T.rds")

all <- readRDS("/hwdata/home/zcshixi7/project/pymt/analyze/Tcell/harmony/T.rds")
FeaturePlot(all,reduction = "umap",features = c("MKI67","STMN1","MCM3"),ncol = 2,raster = F)
DimPlot(all, reduction = "umap",label = TRUE,raster = F)
DimPlot(all, reduction = "umap",label = TRUE,raster = F,group.by = "RNA_snn_res.1.5")
DimPlot(all, reduction = "umap",label = TRUE,raster = F,group.by = "expansion")

table(all$orig.ident,all$Tissue)

DefaultAssay(all) <- "RNA"

FeaturePlot(all,reduction = "umap",features = c("TRDV1","TRDC","TRDV2","TRGV9"),ncol = 2,raster = F)
FeaturePlot(all,reduction = "umap",features = c("nFeature_RNA","nCount_RNA","percent.mt","S.Score","G2M.Score"),ncol = 3,raster = F)


VlnPlot(all,features = c("SELL","CCR7","GZMK","IL7R","GPR183","ITGAE"),group.by = "RNA_snn_res.2",ncol = 1,pt.size = 0)
FeaturePlot(all,reduction = "umap",features = c("CD3E","CD3D","CD4","CD8A","CD3G","TRAC","C1QC","CCR7","TCF7"),ncol =3,raster = F)
FeaturePlot(all,reduction = "umap",features = c("PDCD1","LAG3","HAVCR2","TIGIT","FOXP3","CCR7","TCF7","GZMK","IL7R"),ncol =3,raster = F)


all <- subset(all, idents = c(7, 13), invert = TRUE)

new.cluster.ids <- c("0"="CD8_Trm", 
                     "1"="CD8_Tem", 
                     "2"="CD4_Tn", 
                     "3"="CD8_Teff", 
                     "4"="CD4_Tem", 
                     "5"="CD8_Tex", 
                     "6"="CD4_Tem", 
                     "7"="CD4_Treg", 
                     "8"="CD8_Tex", 
                     "9"="CD4_Tex", 
                     "10"="CD4_Treg", 
                     "11"="CD4_Tex", 
                     "12"="NK", 
                     "13"="gdT", 
                     "14"="NK", 
                     "15"="CD4_Tem", 
                     "16"="CD4_Tn", 
                     "17"="CD4_Treg",
                     "18"="CD4_Tn", 
                     "19"="CD4_Tem", 
                     "20"="Proliferating",
                     "21"="gdT_G9D2", 
                     "22"="CD8_Teff", 
                     "23"="CD8_Tem",
                     "24"="CD4_Treg",
                     "25"="CD8_Tem", 
                     "26"="non-Tcell", 
                     "27"="CD8_Tem",
                     "28"="CD8_Teff", 
                     "29"="CD8_Tn", 
                     "30"="Proliferating",
                     "31"="CD8_Tex", 
                     "32"="non-Tcell", 
                     "33"="CD8_Trm",
                     "34"="Proliferating",
                     "35"="CD4_Tn", 
                     "36"="NK",
                     "37"="Proliferating"
)

# Add cell type to the metadata as a new column
all <- RenameIdents(all, new.cluster.ids)
all$celltype <- all@active.ident

DimPlot(all, reduction = "umap",label = TRUE,raster = F,group.by = "celltype")

saveRDS(all,"T.rds")