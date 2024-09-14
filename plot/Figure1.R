---
title: "PyMT_Figure1"
author: "Chuang Zhao"
date: "2024-09-14"
output: html_document
---

#Merge and Dimensionality Reduction
nonimmune <- readRDS("/home/shixi7/zhaochuang/project/pymt/PyMT/normal-annotion/seurat/nonimmune_annotion.rds")
immune <- readRDS("/home/shixi7/zhaochuang/project/pymt/immune/seurat/decontX_immune/immune_annotion.rds")
all <- merge(x=nonimmune,y=immune)

# Find variable genes for PCA 
all <- FindVariableFeatures(object = all,selection.method = "vst", nfeatures = 2000)
# Scale the data
all.gene <- rownames(all)
all <- ScaleData(all,features = all.gene)

# Perform linear dimensional reduction
all <- RunPCA(object=all,features = VariableFeatures(object = all))
DimPlot(all, reduction = "pca")
ElbowPlot(all,n=30)

# Cluster the cells
all <- FindNeighbors(all, reduction = "pca", dims = 1:30)
all <- FindClusters(all,resolution = 1.5, algorithm = 1)

# Run non-linear dimensional reduction (tSNE and uMAP)
all <- RunTSNE(object=all,dims.use=1:30,do.fast=TRUE,check_duplicates = FALSE)
all <- RunUMAP(all, reduction = "pca", dims = 1:30,min.dist = 1,n.neighbors = 60L)


#Plot Figure 1d
DimPlot(all, reduction = "umap",label = TRUE,group.by="celltype",split.by = "intervene")

all@meta.data$celltype <-
  ordered(all@meta.data$celltype, levels = c("Epithelial","SMC","Fibroblast","Pericyte","Endothelial",
                                             "Tcell","NK","Bcell","Monocyte","Macrophage","Neutrophil","DC"))

use_colors <- c(Epithelial="#D73027",SMC="#F46D43",Fibroblast="#FEB24C",Pericyte="#FEE08B",Endothelial="#FCBBA1",
                Tcell="#6BAED6",NK="#08519C",Bcell="#99D8C9",Monocyte="#7FBC41",Macrophage="#006D2C",Neutrophil="#66C2A4",DC="#00441B")

p1 <- DimPlot(all, reduction = "umap",label = F,
              group.by = "celltype",split.by = "intervene",
              cols = use_colors,pt.size = 1,label.size = 8)
p1
ggsave("nonimmune_UMAP.svg", plot = p1, device = "svg",height = 10,width = 20)

#nonimmune
library(ggalluvial)
nonimmune@meta.data$celltype <-
  ordered(nonimmune@meta.data$celltype, levels = c("Epithelial","SMC","Fibroblast","Pericyte","Endothelial"))
Ratio <- nonimmune@meta.data %>%
  group_by(intervene, celltype) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(relative_freq = n/sum(n))
#select the color
p2 <- ggplot(Ratio, aes(x =intervene, y= relative_freq, fill = celltype,
                        stratum=celltype, alluvium=celltype)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.4, knot.pos=0.5)+
  theme_classic() +
  labs(x='intervene',y = 'Ratio')+
  scale_fill_manual(values = use_colors)
p2
ggsave("nonimmune_Sankey_diagram.svg", plot = p2, device = "svg",height = 10,width = 6)

#immune
immune@meta.data$celltype <-
  ordered(immune@meta.data$celltype, levels = c("Tcell","NK","Bcell","Monocyte","Macrophage","Neutrophil","DC"))
Ratio <- immune@meta.data %>%
  group_by(intervene, celltype) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(relative_freq = n/sum(n))
#select the color
p3 <- ggplot(Ratio, aes(x =intervene, y= relative_freq, fill = celltype,
                        stratum=celltype, alluvium=celltype)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.4, knot.pos=0.5)+
  theme_classic() +
  labs(x='intervene',y = 'Ratio')+
  scale_fill_manual(values = use_colors)
p3
ggsave("immune_Sankey_diagram.svg", plot = p3, device = "svg",height = 10,width = 6)


#Plot Figure 1e
Idents(all) <- "celltype"
# Find markers 
markers <- FindAllMarkers(object = all, 
                          only.pos = TRUE,
                          logfc.threshold= 0.25,
                          min.pct = 0.25)

new <- c("Epithelial","SMC","Fibroblast","Pericyte","Endothelial",
         "Tcell","NK","Bcell","Monocyte","Macrophage","Neutrophil","DC")
new1 <- rev(new)
new1

all@meta.data$celltype <-ordered(all@meta.data$celltype, levels = new1)
markers <- subset(markers,gene%in%c("Epcam","Krt18","PyMTexo",
                                    "Acta2","Tagln","Myl9",
                                    "Col1a1","Dcn","Lum",
                                    "Rgs5","Pdgfrb","Mcam",
                                    "Pecam1","Plvap","Cldn5",
                                    "Cd3g","Cd3e","Cd3d",
                                    "Ncr1","Klrd1","Nkg7",
                                    "Cd79a","Ms4a1","Cd79b",
                                    "Lyz2","Plac8","Cd14",
                                    "Apoe","C1qc","C1qa",
                                    "S100a8","S100a9","Acod1",
                                    "H2-Aa","H2-Eb1","H2-Ab1"
))

row.names(markers) <- 1:nrow(markers)
#Select the row to show
markers2 <- markers[c(select-the-row-to-show),]

Idents(all) <- "celltype"

p <- DotPlot(all,
             features = split(markers2$gene, markers2$cluster),
             cols = c("#ffffff", "#448444")
) +
  RotatedAxis() + # 来自Seurat
  theme(
    panel.border = element_rect(color = "black"),
    panel.spacing = unit(1, "mm"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
  )
p
p$data$feature.groups2 <- factor(p$data$feature.groups,levels = new)
p$data$features.plot <- factor(p$data$features.plot,levels = c("Epcam","Krt18","PyMTexo",
                                                               "Acta2","Tagln","Myl9",
                                                               "Col1a1","Dcn","Lum",
                                                               "Rgs5","Pdgfrb","Mcam",
                                                               "Pecam1","Plvap","Cldn5",
                                                               "Cd3g","Cd3e","Cd3d",
                                                               "Ncr1","Klrd1","Nkg7",
                                                               "Cd79a","Ms4a1","Cd79b",
                                                               "Lyz2","Plac8","Cd14",
                                                               "Apoe","C1qc","C1qa",
                                                               "S100a8","S100a9","Acod1",
                                                               "H2-Aa","H2-Eb1","H2-Ab1"))

library(ggh4x)
library(ggplot2)
library(dplyr)
strip <- strip_themed(background_x = elem_list_rect(fill =c("#D73027","#F46D43","#FEB24C","#FEE08B","#FCBBA1","#6BAED6","#08519C","#99D8C9","#7FBC41","#006D2C","#66C2A4","#00441B") ))

p <- p$data %>% 
  ggplot(aes(x = features.plot,
             y = id)) + 
  geom_point(aes(size = pct.exp, 
                 color = avg.exp.scaled)) + 
  facet_wrap2(~feature.groups2, 
              scales = "free_x", 
              strip = strip, 
              nrow = 1) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 0.5, 
                                   vjust = 0.3, 
                                   color = "black"),
        axis.title = element_blank(),
        strip.background = element_rect(color = "white"),
        axis.text.y = element_blank()) + 
  scale_color_gradient(low = "white",
                       high = "#E31A1C", 
                       name = "avg.exp")
p


df <- data.frame(x = 0, y = levels(all), stringsAsFactors = F )
df$y <- factor(df$y, levels = df$y )

pl <- ggplot(df, aes(x, y, color = factor(y))) +
  geom_point(size = 6, shape = 15,show.legend = F) +
  scale_color_manual(values = rev(c("#D73027","#F46D43","#FEB24C","#FEE08B","#FCBBA1","#6BAED6","#08519C","#99D8C9","#7FBC41","#006D2C","#66C2A4","#00441B"))) +
  theme_classic() +
  scale_x_continuous(expand = c(0,0)) +
  theme(
    plot.margin = margin(r=0),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 9),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
pl

library(cowplot)
p3 <- plot_grid(pl, p, align = "h", axis="bt", rel_widths = c(1.5, 9))
p3

ggsave("dotplot.svg", plot = p3, device = "svg",height = 6,width = 13)















