# ------------------------------------------------------------------
# Integration and Dimensionality Reduction of Myeloid Cells in the PyMT Mouse Model
# ------------------------------------------------------------------
#Convert the Human Gene List
convertHumanGeneList <- function(x){
  require("biomaRt")
  human <- biomaRt::useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
  mouse <- biomaRt::useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
m.s.genes <- convertHumanGeneList(cc.genes$s.genes)
m.g2m.genes <- convertHumanGeneList(cc.genes$g2m.genes)

#Myeloid Cells Reduction
immune <- readRDS("../immune_annotion.rds")
unique(immune$celltype)
Myeloid <- subset(immune,celltype%in%c("Monocyte","Macrophage","Neutrophil","DC"))
rm(immune)

# Find variable genes for PCA 
Myeloid <- FindVariableFeatures(object = Myeloid,selection.method = "vst", nfeatures = 2000)

#Assign a cell cycle score
Myeloid <- CellCycleScoring(Myeloid, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE) 
Myeloid <- RunPCA(Myeloid, features = c(m.s.genes, m.g2m.genes))
DimPlot(Myeloid,reduction = 'pca')

#Remove cycle effects
Myeloid <- ScaleData(Myeloid,
                     vars.to.regress = c("S.Score", "G2M.Score"),
                     features = rownames(Myeloid))

# Perform linear dimensional reduction
Myeloid <- RunPCA(object=Myeloid,features = VariableFeatures(object = Myeloid))
DimPlot(Myeloid, reduction = "pca")
ElbowPlot(Myeloid,n=30)

# Cluster the cells
Myeloid <- FindNeighbors(Myeloid, reduction = "pca", dims = 1:30)
Myeloid <- FindClusters(Myeloid,resolution = 1.5, algorithm = 1)

# Run non-linear dimensional reduction (tSNE and uMAP)
Myeloid <- RunTSNE(object=Myeloid,dims.use=1:30,do.fast=TRUE,check_duplicates = FALSE)
Myeloid <- RunUMAP(Myeloid, reduction = "pca", dims = 1:30)

DimPlot(Myeloid, reduction = "umap",label = TRUE)

# Find markers 
markers <- FindAllMarkers(object = Myeloid, 
                          only.pos = TRUE,
                          logfc.threshold= 0.25,
                          min.pct = 0.25)

#write the markers of each cluster
write.table(markers,file="Myeloid_cycling_dims30_resolution1.5.tsv",sep="\t")

Idents(Myeloid) <- "seurat_clusters"
new.cluster.ids <- c("0"="Macrophage", 
                     "1"="Macrophage", 
                     "2"="Neutrophil", 
                     "3"="Monocyte", 
                     "4"="Neutrophil", 
                     "5"="pre_DC", 
                     "6"="DC2", 
                     "7"="Macrophage", 
                     "8"="Neutrophil", 
                     "9"="Macrophage", 
                     "10"="DC1", 
                     "11"="Macrophage", 
                     "12"="pDC", 
                     "13"="DC3"
)
# Add cell type to the metadata as a new column
Myeloid <- RenameIdents(Myeloid, new.cluster.ids)
Myeloid$subcelltype1 <- Myeloid@active.ident

#Save the RDS
saveRDS(Myeloid,"Myeloid_annotion.rds")


#Annotion the MPS
Myeloid <- readRDS("../Myeloid_annotion.rds")
unique(Myeloid$subcelltype1)
MPS <- subset(Myeloid,subcelltype1%in%c("Macrophage","Monocyte"))
rm(Myeloid)

# Find variable genes for PCA 
MPS <- FindVariableFeatures(object = MPS,selection.method = "vst", nfeatures = 2000)

#Assign a cell cycle score
MPS <- CellCycleScoring(MPS, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE) 
MPS <- RunPCA(MPS, features = c(m.s.genes, m.g2m.genes))
DimPlot(MPS,reduction = 'pca')

#Remove cycle effects
MPS <- ScaleData(MPS,
                 vars.to.regress = c("S.Score", "G2M.Score"),
                 features = rownames(MPS))

# Perform linear dimensional reduction
MPS <- RunPCA(object=MPS,features = VariableFeatures(object = MPS))
DimPlot(MPS, reduction = "pca")
ElbowPlot(MPS,n=30)

# Cluster the cells
MPS <- FindNeighbors(MPS, reduction = "pca", dims = 1:30)
MPS <- FindClusters(MPS,resolution = 0.8, algorithm = 1)

# Run non-linear dimensional reduction (tSNE and uMAP)
MPS <- RunTSNE(object=MPS,dims.use=1:30,do.fast=TRUE,check_duplicates = FALSE)
MPS <- RunUMAP(MPS, reduction = "pca", dims = 1:30)
DimPlot(MPS, reduction = "umap",label = TRUE)

#Annotion
Idents(MPS) <- "seurat_clusters"
new.cluster.ids <- c("0"="Mac_Ifnb1", 
                     "1"="Mac_H2_Ab1", 
                     "2"="Mono_Ccr2", 
                     "3"="Mac_Spp1_Vegfa", 
                     "4"="Mac_Slc40a1", 
                     "5"="Mac_Fcrls"
)
# Add cell type to the metadata as a new column
MPS <- RenameIdents(MPS, new.cluster.ids)
MPS$subcelltype2 <- MPS@active.ident

#Macrophages Reduction
Macro <- subset(MPS,subcelltype2%in%c("Mac_H2_Ab1","Mac_Slc40a1","Mac_Ifnb1","Mac_Spp1_Vegfa","Mac_Fcrls"))

# Find variable genes for PCA
Macro <- FindVariableFeatures(object = Macro,selection.method = "vst", nfeatures = 2000)

#Assign a cell cycle score
Macro <- CellCycleScoring(Macro, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE) 
Macro <- RunPCA(Macro, features = c(m.s.genes, m.g2m.genes))
DimPlot(Macro,reduction = 'pca')

#Remove cycle effects
Macro <- ScaleData(Macro,
                   vars.to.regress = c("S.Score", "G2M.Score"),
                   features = rownames(Macro))

# Perform linear dimensional reduction
Macro <- RunPCA(object=Macro,features = VariableFeatures(object = Macro))
DimPlot(Macro, reduction = "pca")
ElbowPlot(Macro,n=30)

# Cluster the cells
Macro <- FindNeighbors(Macro, reduction = "pca", dims = 1:30)
Macro <- FindClusters(Macro,resolution = 1.2, algorithm = 1)

# Run non-linear dimensional reduction (tSNE and uMAP)
Macro <- RunTSNE(object=Macro,dims.use=1:30,do.fast=TRUE,check_duplicates = FALSE)
Macro <- RunUMAP(Macro, reduction = "pca", dims = 1:30)

DimPlot(Macro, reduction = "umap",label = TRUE)

#Find markers 
markers <- FindAllMarkers(object = Macro, 
                          only.pos = TRUE,
                          logfc.threshold= 0.25,
                          min.pct = 0.25)

#write the markers of each cluster
write.table(markers,file="Macro_cycling_dims30_resolution1.5.tsv",sep="\t")

#Annotion
Idents(Macro) <- "seurat_clusters"
new.cluster.ids <- c("0"="Mac_Ccr2", 
                     "1"="Mac_Ifnb1", 
                     "2"="Mac_Ifnb1", 
                     "3"="Mac_Vegfa", 
                     "4"="Mac_Jun", 
                     "5"="Mac_Gdf15",
                     "6"="Mac_Fcrls"
)
# Add cell type to the metadata as a new column
Macro <- RenameIdents(Macro, new.cluster.ids)
Macro$subcelltype2 <- Macro@active.ident

#Annotion
meta <- MPS@meta.data
meta$cell_name <- rownames(meta)
meta$cluster <- as.character(MPS@active.ident)
ident <- meta$cluster
names(ident) <- rownames(meta)

ident1 <- as.character(Macro@active.ident)
names(ident1) <- names(Macro@active.ident)
for(i in names(ident1)){
  if( i %in% meta$cell_name){
    ident[i] <- ident1[i]}
}
meta$subcelltype2 <- ident
MPS@meta.data <- meta
Idents(MPS) <-'subcelltype2'

#Save the RDS
saveRDS(MPS,"MPS_annotion.rds")


#Neutrophil Reduction
Myeloid <- readRDS("../Myeloid_annotion.rds")
unique(Myeloid$subcelltype1)
Neutrophil <- subset(Myeloid,subcelltype1%in%c("Neutrophil"))

# Find variable genes for PCA 
Neutrophil <- FindVariableFeatures(object = Neutrophil,selection.method = "vst", nfeatures = 2000)

#Assign a cell cycle score
Neutrophil <- CellCycleScoring(Neutrophil, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE) 
Neutrophil <- RunPCA(Neutrophil, features = c(m.s.genes, m.g2m.genes))
DimPlot(Neutrophil,reduction = 'pca')

#Remove cycle effects
Neutrophil <- ScaleData(Neutrophil,
                        vars.to.regress = c("S.Score", "G2M.Score"),
                        features = rownames(Neutrophil))

# Perform linear dimensional reduction
Neutrophil <- RunPCA(object=Neutrophil,features = VariableFeatures(object = Neutrophil))
DimPlot(Neutrophil, reduction = "pca")
ElbowPlot(Neutrophil,n=30)

# Cluster the cells
Neutrophil <- FindNeighbors(Neutrophil, reduction = "pca", dims = 1:30)
Neutrophil <- FindClusters(Neutrophil,resolution = 1.2, algorithm = 1)

# Run non-linear dimensional reduction (tSNE and uMAP)
Neutrophil <- RunTSNE(object=Neutrophil,dims.use=1:30,do.fast=TRUE,check_duplicates = FALSE)
Neutrophil <- RunUMAP(Neutrophil, reduction = "pca", dims = 1:30)

DimPlot(Neutrophil, reduction = "umap",label = TRUE)

#Annotion
Idents(Neutrophil) <- "seurat_clusters"
new.cluster.ids <- c("0"="Neu_Sell", 
                     "1"="Neu_Icam1", 
                     "2"="Neu_Vegfa",
                     "3"="Neu_Cxcl10"
)
# Add cell type to the metadata as a new column
Neutrophil <- RenameIdents(Neutrophil, new.cluster.ids)
Neutrophil$subcelltype2 <- Neutrophil@active.ident

# Find markers 
markers <- FindAllMarkers(object = Neutrophil, 
                          only.pos = TRUE,
                          logfc.threshold= 0.25,
                          min.pct = 0.25)

#write the markers of each cluster
write.table(markers,file="Neutrophil_cycling_dims30_resolution0.8.tsv",sep="\t")

#Save the RDS
saveRDS(Neutrophil,"Neutrophil_annotion.rds")


#Map the MPS to Myeloid
meta <- Myeloid@meta.data
meta$cell_name <- rownames(meta)
meta$cluster <- as.character(Myeloid@active.ident)
ident <- meta$cluster
names(ident) <- rownames(meta)

ident1 <- as.character(MPS@active.ident)
names(ident1) <- names(MPS@active.ident)
for(i in names(ident1)){
  if( i %in% meta$cell_name){
    ident[i] <- ident1[i]}
}
meta$subcelltype2 <- ident
Myeloid@meta.data <- meta
Idents(Myeloid) <-'subcelltype2'

#Map the Neutrophil to Myeloid
meta <- Myeloid@meta.data
meta$cell_name <- rownames(meta)
meta$cluster <- as.character(Myeloid@active.ident)
ident <- meta$cluster
names(ident) <- rownames(meta)

ident1 <- as.character(Neutrophil@active.ident)
names(ident1) <- names(Neutrophil@active.ident)
for(i in names(ident1)){
  if( i %in% meta$cell_name){
    ident[i] <- ident1[i]}
}
meta$subcelltype2 <- ident
Myeloid@meta.data <- meta
Idents(Myeloid) <-'subcelltype2'

Myeloid@meta.data$subcelltype2[Myeloid@meta.data$subcelltype1 == "pDC"] <- "pDC_Siglech"
Myeloid@meta.data$subcelltype2[Myeloid@meta.data$subcelltype1 == "pre_DC"] <- "pre_DC_Flt3"
Myeloid@meta.data$subcelltype2[Myeloid@meta.data$subcelltype1 == "DC1"] <- "DC1_Xcr1"
Myeloid@meta.data$subcelltype2[Myeloid@meta.data$subcelltype1 == "DC2"] <- "DC2_Cd209a"
Myeloid@meta.data$subcelltype2[Myeloid@meta.data$subcelltype1 == "DC3"] <- "DC3_Ccr7"

#Plot and Save the RDS
DimPlot(Myeloid, reduction = "umap",label = TRUE,group.by = "subcelltype2")
saveRDS(Myeloid,"PyMT_Myeloid.rds")
