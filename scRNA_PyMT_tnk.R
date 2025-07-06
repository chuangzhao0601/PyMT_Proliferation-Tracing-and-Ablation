# ------------------------------------------------------------------
# Integration and Dimensionality Reduction of T/NK Cells in the PyMT Mouse Model
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

#Tcell/NK Reduction
immune <- readRDS("/home/shixi7/zhaochuang/project/pymt/immune/seurat/decontX_immune/immune_annotion.rds")
Tcell <- subset(immune,celltype%in%c("Tcell","NK"))
rm(immune)

# Find variable genes for PCA 
Tcell <- FindVariableFeatures(object = Tcell,selection.method = "vst", nfeatures = 2000)

#Assign a cell cycle score
Tcell <- CellCycleScoring(Tcell, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE) 
Tcell <- RunPCA(Tcell, features = c(m.s.genes, m.g2m.genes))
DimPlot(Tcell,reduction = 'pca')

#Remove cycle effects
Tcell <- ScaleData(Tcell,
                   vars.to.regress = c("S.Score", "G2M.Score"),
                   features = rownames(Tcell))

# Perform linear dimensional reduction
Tcell <- RunPCA(object=Tcell,features = VariableFeatures(object = Tcell))
DimPlot(Tcell, reduction = "pca")
ElbowPlot(Tcell,n=30)

# Cluster the cells
Tcell <- FindNeighbors(Tcell, reduction = "pca", dims = 1:30)
Tcell <- FindClusters(Tcell,resolution = 1.5, algorithm = 1)

# Run non-linear dimensional reduction (tSNE and uMAP)
Tcell <- RunTSNE(object=Tcell,dims.use=1:30,do.fast=TRUE,check_duplicates = FALSE)
Tcell <- RunUMAP(Tcell, reduction = "pca", dims = 1:30)

DimPlot(Tcell, reduction = "umap",label = TRUE)

# Find markers 
markers <- FindAllMarkers(object = Tcell, 
                          only.pos = TRUE,
                          logfc.threshold= 0.25,
                          min.pct = 0.25)

#write the markers of each cluster
write.table(markers,file="Tcell_cycling_dims30_resolution1.5.tsv",sep="\t")

Idents(Tcell) <- "seurat_clusters"
new.cluster.ids <- c("0"="CD8", 
                     "1"="CD8", 
                     "2"="CD8", 
                     "3"="CD4", 
                     "4"="CD8", 
                     "5"="CD8", 
                     "6"="NK", 
                     "7"="Treg", 
                     "8"="γδT", 
                     "9"="NK", 
                     "10"="Th17", 
                     "11"="CD8", 
                     "12"="CD8", 
                     "13"="NK", 
                     "14"="CD8"
)
# Add cell type to the metadata as a new column
Tcell <- RenameIdents(Tcell, new.cluster.ids)
Tcell$subcelltype1 <- Tcell@active.ident

#saveRDS(Tcell,"annotion_Tcell.rds")


#CD8T/NK Reduction
Tcell <- readRDS("/home/shixi7/zhaochuang/project/pymt/immune/seurat/decontX_immune/T_NK/annotion_Tcell.rds")
DimPlot(Tcell, reduction = "umap",label = TRUE)
unique(Tcell$subcelltype1)
CD8 <- subset(Tcell,subcelltype1%in%c("CD8","NK"))

# Find variable genes for PCA 
CD8 <- FindVariableFeatures(object = CD8,selection.method = "vst", nfeatures = 2000)

#Assign a cell cycle score
CD8 <- CellCycleScoring(CD8, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE) 
CD8 <- RunPCA(CD8, features = c(m.s.genes, m.g2m.genes))
DimPlot(CD8,reduction = 'pca')

#Remove cycle effects
CD8 <- ScaleData(CD8,vars.to.regress = c("S.Score", "G2M.Score"),features = rownames(CD8))

# Perform linear dimensional reduction
CD8 <- RunPCA(object=CD8,features = VariableFeatures(object = CD8))
DimPlot(CD8, reduction = "pca")
ElbowPlot(CD8,n=30)

# Cluster the cells
CD8 <- FindNeighbors(CD8, reduction = "pca", dims = 1:30)
CD8 <- FindClusters(CD8,resolution = 1.5, algorithm = 1)

# Run non-linear dimensional reduction (tSNE and uMAP)
CD8 <- RunTSNE(object=CD8,dims.use=1:30,do.fast=TRUE,check_duplicates = FALSE)
CD8 <- RunUMAP(CD8, reduction = "pca", dims = 1:30)

DimPlot(CD8, reduction = "umap",label = TRUE)

# Find markers 
markers <- FindAllMarkers(object = CD8, 
                          only.pos = TRUE,
                          logfc.threshold= 0.25,
                          min.pct = 0.25)

#write the markers of each cluster
write.table(markers,file="CD8_cycling_dims30_resolution1.5.tsv",sep="\t")

Idents(CD8) <- "seurat_clusters"
new.cluster.ids <- c("0"="CD8_Pdcd1", 
                     "1"="NKT_Gzmb", 
                     "2"="CD8_Gzmk", 
                     "3"="CD8_Tcf7", 
                     "4"="NKT_Tcf7",
                     "5"="CD8_Ifit3", 
                     "6"="NKT_Gzmb", 
                     "7"="NK_Ccr2", 
                     "8"="NK_Gzmc", 
                     "9"="CD8_Ifit3",
                     "10"="CD8_Gzmk", 
                     "11"="NK_Sell", 
                     "12"="CD8_Ccr7", 
                     "13"="CD8_Hspa1a"
)
# Add cell type to the metadata as a new column
CD8 <- RenameIdents(CD8, new.cluster.ids)
CD8$subcelltype2 <- CD8@active.ident

#Save the RDS
saveRDS(CD8,"CD8_cycling.rds")

#Rename the CD4T
CD4 <- subset(Tcell,subcelltype1%in%c("CD4","Th17","Treg"))
DimPlot(CD4, reduction = "umap",label = TRUE)

Idents(CD4) <- "subcelltype1"
new.cluster.ids <- c("CD4"="CD4_Tcf7", 
                     "Th17"="Th17_Il17a", 
                     "Treg"="Treg_Foxp3"
)
# Add cell type to the metadata as a new column
CD4 <- RenameIdents(CD4, new.cluster.ids)
CD4$subcelltype2 <- CD4@active.ident

#Annotion CD8
meta <- Tcell@meta.data
meta$cell_name <- rownames(meta)
meta$cluster <- as.character(Tcell@active.ident)
ident <- meta$cluster
names(ident) <- rownames(meta)
ident1 <- as.character(CD8@active.ident)
names(ident1) <- names(CD8@active.ident)
for(i in names(ident1)){
  if( i %in% meta$cell_name){
    ident[i] <- ident1[i]}
}
meta$subcelltype2 <- ident
Tcell@meta.data <- meta
Idents(Tcell) <-'subcelltype2'

#Annotion CD4
meta <- Tcell@meta.data
meta$cell_name <- rownames(meta)
meta$cluster <- as.character(Tcell@active.ident)
ident <- meta$cluster
names(ident) <- rownames(meta)
ident1 <- as.character(CD4@active.ident)
names(ident1) <- names(CD4@active.ident)
for(i in names(ident1)){
  if( i %in% meta$cell_name){
    ident[i] <- ident1[i]}
}
meta$subcelltype2 <- ident
Tcell@meta.data <- meta
Idents(Tcell) <-'subcelltype2'

DimPlot(Tcell, reduction = "umap",label = TRUE)

#Save the RDS
saveRDS(Tcell,"Tcell_final_annotion.rds")