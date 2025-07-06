# ------------------------------------------------------------------
# cellchat analysis on tumor and immune cells (mouse data)
# ------------------------------------------------------------------
library(CellChat)
library(Seurat)
library(future)
plan("multicore", workers = 4)
setwd("/home/shixi7/zhaochuang/project/pymt/cellchat_cellphone/cellchat/")

# ---------------------------------------------------------------------------
# 1. Load Seurat object and set cell type
# ---------------------------------------------------------------------------
all <- readRDS("/home/shixi7/zhaochuang/project/pymt/PyMT/normal-annotion/seurat/Epithelial/all.rds")

# Define cell type order
ctype_levels <- c(
  "Epithelial", "SMC", "Fibroblast", "Pericyte", "Endothelial",
  "Tcell", "NK", "Bcell", "Monocyte", "Macrophage", "Neutrophil", "DC"
)
all$celltype <- factor(all$celltype, levels = ctype_levels)
Idents(all) <- "celltype"

# Optional: visualize UMAP colored by cell type
DimPlot(all, reduction = "umap", group.by = "celltype", label = TRUE)

# ---------------------------------------------------------------------------
# 2. Create CellChat object
# ---------------------------------------------------------------------------
all$cell_annotations <- Idents(all)
cellchat <- createCellChat(object = all)
groupSize <- as.numeric(table(cellchat@idents))

# ---------------------------------------------------------------------------
# 3. Load ligand-receptor database
# ---------------------------------------------------------------------------
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
cellchat@DB <- subsetDB(CellChatDB)

# ---------------------------------------------------------------------------
# 4. Preprocessing: gene filtering & projection
# ---------------------------------------------------------------------------
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)

# ---------------------------------------------------------------------------
# 5. Communication probability inference
# ---------------------------------------------------------------------------
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)

# ---------------------------------------------------------------------------
# 6. Inference at pathway level
# ---------------------------------------------------------------------------
cellchat <- computeCommunProbPathway(cellchat)

# ---------------------------------------------------------------------------
# 7. Visualization: example for "CXCL" pathway
# ---------------------------------------------------------------------------
pathways.show <- "CXCL"
group.colors <- c(
  Epithelial = "#D73027", SMC = "#F46D43", Fibroblast = "#FEB24C", Pericyte = "#FEE08B",
  Endothelial = "#FCBBA1", Tcell = "#6BAED6", NK = "#08519C", Bcell = "#99D8C9",
  Monocyte = "#7FBC41", Macrophage = "#006D2C", Neutrophil = "#66C2A4", DC = "#00441B"
)

# Circle plot
svg("CXCL_circle.svg", width = 8, height = 8)
netVisual_aggregate(cellchat, signaling = pathways.show,
                    layout = "circle", arrow.size = 0.8,
                    color.use = group.colors)
dev.off()

# Chord plot
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

# ---------------------------------------------------------------------------
# 8. Visualize individual ligand-receptor pair
# ---------------------------------------------------------------------------
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[2, ]  # example pair

cellchat@idents <- factor(as.character(cellchat@idents), levels = ctype_levels)
cellchat@meta$ident <- cellchat@idents
vertex.receiver <- seq_len(length(ctype_levels))  # mark all as receiver

p1 <- netVisual_individual(cellchat,
                           signaling = pathways.show,
                           pairLR.use = LR.show,
                           vertex.receiver = vertex.receiver,
                           arrow.size = 0.8,
                           color.use = group.colors)
ggsave("CXCL2-CXCR2.svg", p1, height = 5, width = 5)

# ---------------------------------------------------------------------------
# 9. Visualize gene expression
# ---------------------------------------------------------------------------
p2 <- plotGeneExpression(cellchat, signaling = "CXCL")
ggsave("CXCL_expression_allGenes.pdf", p2, width = 10, height = 6)

p3 <- plotGeneExpression(cellchat, signaling = "CXCL",
                         features = c("Cxcl2", "Cxcr2"),
                         color.use = group.colors)
ggsave("CXCL_expression_selected.svg", p3, width = 11, height = 6)

# ---------------------------------------------------------------------------
# 10. Save CellChat object
# ---------------------------------------------------------------------------
saveRDS(cellchat, file = "cellchat_CXCL.rds")
