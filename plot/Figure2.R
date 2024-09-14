---
title: "PyMT_Figure2"
author: "Chuang Zhao"
date: "2024-09-14"
output: html_document
---

#Plot Figure 2a
Epithelial <- readRDS("/home/shixi7/zhaochuang/project/pymt/PyMT/normal-annotion/seurat/Epithelial/Epithelial_annotion.rds")
##UMAP
Epithelial@meta.data$subcelltype1 <-
  ordered(Epithelial@meta.data$subcelltype1, levels = c("Cancer_cells_c1_Spp1","Cancer_cells_c2_Tspan8",
                                                        "Cancer_cells_c3_Hspa1a","Cancer_cells_c4_Cxcl1",
                                                        "Cancer_cells_c5_Lalba","Cancer_cells_c6_Pde4d",
                                                        "Cancer_cells_c7_Myl9","Cancer_cells_c8_Ntrk2",
                                                        "Cancer_cells_c9_Top2a","Cancer_cells_c10_Ly6a",
                                                        "Cancer_cells_c11_Ifit3","Cancer_cells_c12_Aldh1a3",
                                                        "Cancer_cells_c13_Prlr","Cancer_cells_c14_Csn1s1"))

use_colors <- c(Cancer_cells_c1_Spp1="#2171B5",Cancer_cells_c2_Tspan8="#F46D43",
                Cancer_cells_c3_Hspa1a="#D73027",Cancer_cells_c4_Cxcl1="#FC9272",
                Cancer_cells_c5_Lalba="#B3DE69",Cancer_cells_c6_Pde4d="#08519C",
                Cancer_cells_c7_Myl9="#FDAE61",Cancer_cells_c8_Ntrk2="#66C2A4",
                Cancer_cells_c9_Top2a="#41AE76",Cancer_cells_c10_Ly6a="#BEBADA",
                Cancer_cells_c11_Ifit3="#FCCDE5",Cancer_cells_c12_Aldh1a3="#9E0142",
                Cancer_cells_c13_Prlr="#5E4FA2",Cancer_cells_c14_Csn1s1="#FEE08B")

p <- DimPlot(Epithelial, reduction = 'umap',group.by = "subcelltype1",split.by = 'intervene', label = F, cols = use_colors,pt.size = 1)
p
ggsave("Epithelial_UMAP.svg", plot = p, device = "svg",height = 5,width = 10)

##Sankey
library(ggalluvial)
Ratio <- Epithelial@meta.data %>%
  group_by(intervene, subcelltype1) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(relative_freq = n/sum(n))
#select the color

p <- ggplot(Ratio, aes(x =intervene, y= relative_freq, fill = subcelltype1,
                        stratum=subcelltype1, alluvium=subcelltype1)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.4, knot.pos=0.5)+
  theme_classic() +
  labs(x='intervene',y = 'Ratio')+
  scale_fill_manual(values = use_colors)
p
ggsave("Epithelial_Sankey_diagram.svg", plot = p, device = "svg",height = 10,width = 10)


#Plot Figure 2b
#Prepare the matrix for scVelo
setwd("/home/shixi7/zhaochuang/project/pymt/PyMT/normal-annotion/seurat/scVelo/")
#提取坐标信息
write.csv(Embeddings(Epithelial, reduction = "umap"), file = "cell_embeddings.csv")
# 获取每个细胞的barcode
write.csv(Cells(Epithelial), file = "cellID_obs.csv", row.names = FALSE)
# 提取每个细胞的cluster信息
write.csv(Epithelial@meta.data[, 'seurat_clusters', drop = FALSE], file = "cell_clusters.csv")
# 提取每个细胞的celltype信息
write.csv(Epithelial@meta.data[, 'subcelltype1', drop = FALSE], file = "cell_celltype.csv")
library(scales)
# 获取celltype的颜色信息
hue_pal()(length(levels(Epithelial$subcelltype1)))


#Plot Figure 2c
Tcell <- readRDS("/home/shixi7/zhaochuang/project/pymt/immune/seurat/decontX_immune/T_NK/Tcell_final_annotion.rds")
##Umap
Tcell@meta.data$subcelltype2 <-
  ordered(Tcell@meta.data$subcelltype2, levels = c("CD4_Tcf7","Th17_Il17a","Treg_Foxp3",
                                                   "CD8_Ccr7", "CD8_Tcf7", "CD8_Ifit3", "CD8_Gzmk", "CD8_Pdcd1", "CD8_Hspa1a", 
                                                   "NK_Sell","NK_Ccr2","NK_Gzmc",
                                                   "NKT_Tcf7","NKT_Gzmb","γδT"))

use_colors <- c(CD4_Tcf7="#D73027",Th17_Il17a="#F46D43",Treg_Foxp3="#FC9272",
                CD8_Ccr7="#FDAE61",CD8_Tcf7="#E0F3F8",CD8_Ifit3="#ABD9E9",CD8_Gzmk="#6BAED6",CD8_Pdcd1="#08519C",CD8_Hspa1a="#5E4FA2",
                NK_Sell="#66C2A4",NK_Ccr2="#B3DE69",NK_Gzmc="#41AE76",
                NKT_Tcf7="#FEE08B",NKT_Gzmb="#9E0142",
                γδT="black")

p <- DimPlot(Tcell, reduction = 'umap',group.by = "subcelltype2",split.by = 'intervene', label = F, cols = use_colors,pt.size = 1)
p
ggsave("Tcell_UMAP.svg", plot = p, device = "svg",height = 5,width = 10)

##Sankey
library(ggalluvial)
Ratio <- Tcell@meta.data %>%
  group_by(intervene, subcelltype2) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(relative_freq = n/sum(n))
#select the color

p <- ggplot(Ratio, aes(x =intervene, y= relative_freq, fill = subcelltype2,
                        stratum=subcelltype2, alluvium=subcelltype2)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.4, knot.pos=0.5)+
  theme_classic() +
  labs(x='intervene',y = 'Ratio')+
  scale_fill_manual(values = use_colors)
p
ggsave("Tcell_Sankey_diagram.svg", plot = p, device = "svg",height = 10,width = 10)


#Plot Figure 2d
γδT <- subset(Tcell,subcelltype2%in%c("γδT"))
##Anti-tumor
p <- VlnPlot(object = γδT, 
              features = c("Gzma","Gzmb","Prf1","Il2","Tnf","Ifng","Ccr5","Fasl","Tnfsf10","Klrk1","Tnfsf9"),
              group.by="subcelltype2", ncol = 4,pt.size = 0,same.y.lims = TRUE)

p
ggsave("VlnPlot_anti.svg", plot = p, device = "svg",height = 10,width = 10)

##Pro-tumor
p <- VlnPlot(object = γδT, 
              features = c("Il17a","Il4","Il10","Lgals1",
                           "Lgals9","Cd274","Vegfa","Csf2"),
              group.by="subcelltype2", ncol = 4,pt.size = 0,same.y.lims = TRUE)

p
ggsave("VlnPlot_pro.svg", plot = p, device = "svg",height = 10,width = 10)


#Plot Figure 2e
Myeloid <- readRDS("/home/shixi7/zhaochuang/project/pymt/immune/seurat/decontX_immune/Myeloid/Myeloid_final_annotion.rds")
#Umap
Myeloid@meta.data$subcelltype2 <-
  ordered(Myeloid@meta.data$subcelltype2, levels = c("Neu_Sell","Neu_Icam1","Neu_Hilpda_Vegfa","Neu_Cxcl10",
                                                     "DC1_Xcr1", "DC2_Cd209a", "DC3_Ccr7", "pre_DC_Flt3", "pDC_Siglech",
                                                     "Mono_Ccr2","Mac_Spp1_Vegfa","Mac_Ccr2",
                                                     "Mac_Ifnb1","Mac_Jun","Mac_Fcrls","Mac_Gdf15"))
use_colors <- c(Neu_Sell="#A50026",Neu_Icam1="#D73027",Neu_Hilpda_Vegfa="#F46D43",Neu_Cxcl10="#FC4E2A",
                DC1_Xcr1="#FDAE61",DC2_Cd209a="#FEE08B",DC3_Ccr7="#D9EF8B",pre_DC_Flt3="#A6D96A",pDC_Siglech="#66BD63",
                Mono_Ccr2="#006837",Mac_Spp1_Vegfa="#9ECAE1",Mac_Ccr2="#4292C6",
                Mac_Ifnb1="#08306B",Mac_Jun="#807DBA",Mac_Fcrls="#313695",Mac_Gdf15="#88419D")

p <- DimPlot(Myeloid, reduction = 'umap',group.by = "subcelltype2",split.by = 'intervene', label = F, cols = use_colors,pt.size = 1)
p
ggsave("Myeloid_UMAP.svg", plot = p, device = "svg",height = 5,width = 10)

##Sankey
library(ggalluvial)
Ratio <- Myeloid@meta.data %>%
  group_by(intervene, subcelltype2) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(relative_freq = n/sum(n))
p <- ggplot(Ratio, aes(x =intervene, y= relative_freq, fill = subcelltype2,
                        stratum=subcelltype2, alluvium=subcelltype2)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.4, knot.pos=0.5)+
  theme_classic() +
  labs(x='intervene',y = 'Ratio')+
  scale_fill_manual(values = use_colors)
p
ggsave("Myeloid_Sankey_diagram.svg", plot = p, device = "svg",height = 10,width = 10)


#Plot Figure 2e
MPS <- readRDS("./Human_bc2021_MPS.rds")
library(ggalluvial)
Ratio <- MPS@meta.data %>%
  group_by(BC_type_expansion, subcelltype1) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(relative_freq = n/sum(n))
#select the color
use_colors <- c("#6BAED6","#D73027")
p <- ggplot(Ratio, aes(x =BC_type_expansion, y= relative_freq, fill = subcelltype1,
                        stratum=subcelltype1, alluvium=subcelltype1)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.4, knot.pos=0.5)+
  theme_classic() +
  labs(x='intervene',y = 'Ratio')+
  scale_fill_manual(values = use_colors)
p
ggsave("BC_type_expansion_Sankey_diagram.svg", plot = p, device = "svg",height = 10,width = 10)
















