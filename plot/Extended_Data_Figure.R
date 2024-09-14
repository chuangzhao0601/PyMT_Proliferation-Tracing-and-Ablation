#Extended Data Figure
---
title: "PyMT_Extended_Data_Figure"
author: "Chuang Zhao"
date: "2024-09-14"
output: html_document
---

#Plot Extended Data Figure 3a
Epithelial <- readRDS("/home/shixi7/zhaochuang/project/pymt/PyMT/normal-annotion/seurat/Epithelial/Epithelial_annotion.rds")

# Find markers 
Idents(Epithelial) <- "subcelltype1"
markers <- FindAllMarkers(object = Epithelial, 
                          only.pos = TRUE,
                          logfc.threshold= 0.25,
                          min.pct = 0.25)
#Other the subcelltype
new <- c("Cancer_cells_c1_Spp1","Cancer_cells_c2_Tspan8","Cancer_cells_c3_Hspa1a","Cancer_cells_c4_Cxcl1",
         "Cancer_cells_c5_Lalba","Cancer_cells_c6_Pde4d","Cancer_cells_c7_Myl9","Cancer_cells_c8_Ntrk2",
         "Cancer_cells_c9_Top2a","Cancer_cells_c10_Ly6a","Cancer_cells_c11_Ifit3","Cancer_cells_c12_Aldh1a3",
         "Cancer_cells_c13_Prlr","Cancer_cells_c14_Csn1s1")
new1 <- rev(new)

Epithelial@meta.data$subcelltype1 <-ordered(Epithelial@meta.data$subcelltype1, levels = new1)
markers <- subset(markers,gene%in%c("Spp1","Chil1","Tspan8","Aldoa","Hspa1a","Hspa1b","Cxcl1","Cxcl2",
                                    "Lalba","Rbp7","Pde4d","Pde7b","Myl9","Igfbp2","Ntrk2","Gjb2",
                                    "Top2a","Pclaf","Ly6a","Ly6d","Ifit3","Isg15","Aldh1a3","Vegfa",
                                    "Prlr","Large1","Csn1s1","Hepacam2"))

row.names(markers) <- 1:nrow(markers) 

markers2 <- markers[c(1,2,3,4,10,11,14,15,18,19,20,21,22,23,28,
                      29,30,31,33,34,39,40,46,47,50,51,53,54),]

p <- DotPlot(Epithelial,
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
p$data$features.plot <- factor(p$data$features.plot,levels = c("Spp1","Chil1",
                                                               "Tspan8","Aldoa",
                                                               "Hspa1a","Hspa1b",
                                                               "Cxcl1","Cxcl2",
                                                               "Lalba","Rbp7",
                                                               "Pde4d","Pde7b",
                                                               "Myl9","Igfbp2",
                                                               "Ntrk2","Gjb2",
                                                               "Top2a","Pclaf",
                                                               "Ly6a","Ly6d",
                                                               "Ifit3","Isg15",
                                                               "Aldh1a3","Vegfa",
                                                               "Prlr","Large1",
                                                               "Csn1s1","Hepacam2"))

library(ggh4x)
library(ggplot2)
library(dplyr)
strip <- strip_themed(background_x = elem_list_rect(fill =c("#2171B5","#F46D43","#D73027","#FC9272","#B3DE69","#08519C","#FDAE61","#66C2A4","#41AE76","#BEBADA","#FCCDE5","#9E0142","#5E4FA2","#FEE08B") ))
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

df <- data.frame(x = 0, y = levels(Epithelial), stringsAsFactors = F )
df$y <- factor(df$y, levels = df$y )

pl <- ggplot(df, aes(x, y, color = factor(y))) +
  geom_point(size = 6, shape = 15,show.legend = F) +
  scale_color_manual(values = rev(c("#2171B5","#F46D43","#D73027","#FC9272","#B3DE69","#08519C","#FDAE61","#66C2A4","#41AE76","#BEBADA","#FCCDE5","#9E0142","#5E4FA2","#FEE08B"))) +
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
p <- plot_grid(pl, p, align = "h", axis="bt", rel_widths = c(1.5, 9))
p

ggsave("dotplot.svg", plot = p, device = "svg",height = 6,width = 13)


#Plot Extended Data Figure 3b
Epithelial <- readRDS("/home/shixi7/zhaochuang/project/pymt/PyMT/normal-annotion/seurat/Epithelial/Epithelial_decontx.rds")
unique(Epithelial$subcelltype1)
immune <- readRDS("/home/shixi7/zhaochuang/project/pymt/immune/seurat/decontX_immune/immune_annotion.rds")
unique(immune$celltype)
#label the intervene
immune@meta.data$subcelltype1 <- "immune"

all <- merge(Epithelial,immune)
all@meta.data$subcelltype1 <-
  ordered(all@meta.data$subcelltype1, levels = c("Cancer_cells_c1_Spp1","Cancer_cells_c2_Tspan8","Cancer_cells_c3_Hspa1a","Cancer_cells_c4_Cxcl1",
                                                 "Cancer_cells_c5_Lalba","Cancer_cells_c6_Pde4d","Cancer_cells_c7_Myl9","Cancer_cells_c8_Ntrk2",
                                                 "Cancer_cells_c9_Top2a","Cancer_cells_c10_Ly6a","Cancer_cells_c11_Ifit3","Cancer_cells_c12_Aldh1a3",
                                                 "Cancer_cells_c13_Prlr","Cancer_cells_c14_Csn1s1","immune"))

use_colors <- c(Cancer_cells_c1_Spp1="#2171B5",Cancer_cells_c2_Tspan8="#F46D43",
                Cancer_cells_c3_Hspa1a="#D73027",Cancer_cells_c4_Cxcl1="#FC9272",
                Cancer_cells_c5_Lalba="#B3DE69",Cancer_cells_c6_Pde4d="#08519C",
                Cancer_cells_c7_Myl9="#FDAE61",Cancer_cells_c8_Ntrk2="#66C2A4",
                Cancer_cells_c9_Top2a="#41AE76",Cancer_cells_c10_Ly6a="#BEBADA",
                Cancer_cells_c11_Ifit3="#FCCDE5",Cancer_cells_c12_Aldh1a3="#9E0142",
                Cancer_cells_c13_Prlr="#5E4FA2",Cancer_cells_c14_Csn1s1="#FEE08B",immune="black")
library(Seurat)
library(ggplot2)
p <- VlnPlot(object = all, 
              features = c("PyMTexo","Epcam"),
              group.by="subcelltype1",pt.size = 0,cols =use_colors,ncol = 1)

p
ggsave("vlnplot.svg", plot = p, device = "svg",height = 13,width = 6)


#Plot Extended Data Figure 4a
Myeloid <- readRDS("/home/shixi7/zhaochuang/project/pymt/immune/seurat/decontX_immune/Myeloid/Myeloid_final_annotion.rds")
Idents(Myeloid) <- "subcelltype2"
# Find markers 
markers <- FindAllMarkers(object = Myeloid, 
                          only.pos = TRUE,
                          logfc.threshold= 0.25,
                          min.pct = 0.25)
#write the markers of each cluster
new <- c("Neu_Sell","Neu_Icam1","Neu_Hilpda_Vegfa","Neu_Cxcl10","DC1_Xcr1", "DC2_Cd209a", "DC3_Ccr7", "pre_DC_Flt3", 
         "pDC_Siglech","Mono_Ccr2","Mac_Spp1_Vegfa","Mac_Ccr2","Mac_Ifnb1","Mac_Jun","Mac_Fcrls","Mac_Gdf15")
new1 <- rev(new)

Myeloid@meta.data$subcelltype2 <-ordered(Myeloid@meta.data$subcelltype2, levels = new1)
markers <- subset(markers,gene%in%c("Sell","Mmp9","Icam1","Cxcr2","Hilpda","Vegfa","Cxcl10","Isg20",
                                    "Xcr1","Btla","Cd209a","Mgl2","Ccr7","Fscn1","Flt3","Tbc1d4",
                                    "Siglech","Ly6d","Ccr2","Id3","Spp1","Mmp12","Tgfbi","Ocstamp",
                                    "Ifnb1","Mmp13","Jun","Hspa1a","Fcrls","Vcam1","Gdf15","Slc40a1"))

row.names(markers) <- 1:nrow(markers) 
markers2 <- markers[c(2,3,5,7,8,10,13,15,18,19,23,24,27,28,31,32,34,
                      35,38,40,42,46,50,51,53,55,59,60,61,62,65,66),]

Idents(Myeloid) <- "subcelltype2"

p <- DotPlot(Myeloid,
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
p$data$features.plot <- factor(p$data$features.plot,levels = c("Sell","Mmp9",
                                                               "Icam1","Cxcr2",
                                                               "Hilpda","Vegfa",
                                                               "Cxcl10","Isg20",
                                                               "Xcr1","Btla",
                                                               "Cd209a","Mgl2",
                                                               "Ccr7","Fscn1",
                                                               "Flt3","Tbc1d4",
                                                               "Siglech","Ly6d",
                                                               "Ccr2","Id3",
                                                               "Spp1","Mmp12",
                                                               "Tgfbi","Ocstamp",
                                                               "Ifnb1","Mmp13",
                                                               "Jun","Hspa1a",
                                                               "Fcrls","Vcam1",
                                                               "Gdf15","Slc40a1"))


library(ggh4x)
strip <- strip_themed(background_x = elem_list_rect(fill =c("#A50026","#D73027","#F46D43","#FC4E2A","#FDAE61","#FEE08B","#D9EF8B","#A6D96A","#66BD63","#006837","#9ECAE1","#4292C6","#08306B","#807DBA","#313695","#88419D") ))

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

df <- data.frame(x = 0, y = levels(Myeloid), stringsAsFactors = F )
df$y <- factor(df$y, levels = df$y )

pl <- ggplot(df, aes(x, y, color = factor(y))) +
  geom_point(size = 6, shape = 15,show.legend = F) +
  scale_color_manual(values = rev(c("#A50026","#D73027","#F46D43","#FC4E2A","#FDAE61","#FEE08B","#D9EF8B","#A6D96A","#66BD63","#006837","#9ECAE1","#4292C6","#08306B","#807DBA","#313695","#88419D"
  ))) +
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
p <- plot_grid(pl, p, align = "h", axis="bt", rel_widths = c(1.5, 9))
p

ggsave("dotplot.svg", plot = p, device = "svg",height = 6,width = 13)


#Plot Extended Data Figure 4b
#Load the Geneset
wlhgeneset <- readxl::read_excel("/home/shixi7/zhaochuang/reference/geneset/wlh.xlsx")
library(msigdbr)
gene_symbol
kegg_list = split(wlhgeneset$gene_symbol, wlhgeneset$gs_name)
kegg_list[1:2]

#Extract the expression matrix
expr <- AverageExpression(Tcell, assays = "decontXcounts", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]
expr <- as.matrix(expr)
head(expr)

library(GSVA)
#GSVA
kegg_gsva <- gsva(expr, kegg_list, kcdf="Gaussian",method = "gsva",parallel.sz=12,verbose=FALSE)
kegg_gsvas <- data.frame(Genesets=rownames(kegg_gsva), kegg_gsva, check.names = F)
head(kegg_gsvas)

#Order the Score and subcelltype1
rownames(kegg_gsva)
desired_order <- c("Naive","Activation_Effector_function","Exhaustion",
                   "TCR_Signaling","Cytotoxicity","Cytokine_Cytokine_receptor","Chemokine_Chemokine_receptor","Senescence",
                   "Anergy","NFKB_Signaling","Stress_response","MAPK_Signaling","Adhesion","IFN_Response",
                   "Oxidative_phosphorylation","Glycolysis","Fatty_acid_metabolism",
                   "Pro_apoptosis","Anti_apoptosis")
new_kegg_gsva <- kegg_gsva[match(desired_order, rownames(kegg_gsva)),match(desired_colorder, colnames(kegg_gsva))]
desired_colorder <- c("CD4_Tcf7","Th17_Il17a","Treg_Foxp3",
                      "CD8_Ccr7", "CD8_Tcf7", "CD8_Ifit3", "CD8_Gzmk", "CD8_Pdcd1", "CD8_Hspa1a", 
                      "NK_Sell","NK_Ccr2","NK_Gzmc",
                      "NKT_Tcf7","NKT_Gzmb","γδT")

p <- pheatmap::pheatmap(new_kegg_gsva,show_colnames = T,cluster_rows=F,cluster_cols=F,
                        scale = "row",angle_col = "45",
                        color = colorRampPalette(c("#6BAED6","#ABD9E9","#E0F3F8","#FFFFBF","#FEE090","#FDAE61","#F46D43"))(50))
p
ggsave(p,file="GSVA_heatmap.svg",width = 8,height = 7)


#Plot Extended Data Figure 4c
#Load the RDS and Genes
Tcell <- readRDS("/home/shixi7/zhaochuang/project/pymt/immune/seurat/decontX_immune/T_NK/Tcell_final_annotion.rds")
final.markers <- c("Pdcd1","Tox","Tigit","Havcr2","Lag3","Vsir","Btla","Ctla4","Cd28","Cd274","Cd226","Cd96","Klrb1")
genes_to_check <- final.markers
library(stringr)  

th <- theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
p_all_markers <- DotPlot(Tcell, features = genes_to_check,
                         assay='RNA' ,group.by = 'subcelltype2' )  + coord_flip()+th

data<-p_all_markers$data

colnames(data)
colnames(data)<-c("AverageExpression_unscaled","Precent_Expressed","Features","Moduletype","Average_Expression")

p <- ggplot(data,aes(Moduletype,Features,size = Precent_Expressed))+
  geom_point(shape=21,aes(fill= Average_Expression),position =position_dodge(0))+
  theme_minimal()+xlab(NULL)+ylab(NULL) +
  scale_size_continuous(range=c(1,10))+theme_bw()+
  scale_fill_stepsn(colours=colorRampPalette(c("#6BAED6","white","#D73027"))(50))+
  theme(legend.position = "right",legend.box = "vertical",
        legend.margin=margin(t= 0, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_text(color="black",size=16,angle = 90, 
                                    vjust = 0.5, hjust=0.5),
        axis.text.y  = element_text(color="black",size=12),
        legend.text = element_text(size =12,color="black"),
        legend.title = element_text(size =12,color="black"),
        axis.title.y=element_text(vjust=1,  
                                  size=16)
  )+labs(x=" ",y = "Features")
p

ggsave("Therapy_target.svg", plot = p, device = "svg",height = 8,width = 8)


#Plot Extended Data Figure 5a
Myeloid <- readRDS("/home/shixi7/zhaochuang/project/pymt/immune/seurat/decontX_immune/Myeloid/Myeloid_final_annotion.rds")
Idents(Myeloid) <- "subcelltype2"
# Find markers 
markers <- FindAllMarkers(object = Myeloid, 
                          only.pos = TRUE,
                          logfc.threshold= 0.25,
                          min.pct = 0.25)
#Order the subcelltype2
new <- c("Neu_Sell","Neu_Icam1","Neu_Hilpda_Vegfa","Neu_Cxcl10",
         "DC1_Xcr1", "DC2_Cd209a", "DC3_Ccr7", "pre_DC_Flt3", "pDC_Siglech",
         "Mono_Ccr2","Mac_Spp1_Vegfa","Mac_Ccr2",
         "Mac_Ifnb1","Mac_Jun","Mac_Fcrls","Mac_Gdf15")
new1 <- rev(new)

Myeloid@meta.data$subcelltype2 <-ordered(Myeloid@meta.data$subcelltype2, levels = new1)
markers <- subset(markers,gene%in%c("Sell","Mmp9","Icam1","Cxcr2","Hilpda","Vegfa","Cxcl10","Isg20",
                                    "Xcr1","Btla","Cd209a","Mgl2","Ccr7","Fscn1","Flt3","Tbc1d4",
                                    "Siglech","Ly6d","Ccr2","Id3","Spp1","Mmp12","Tgfbi","Ocstamp",
                                    "Ifnb1","Mmp13","Jun","Hspa1a","Fcrls","Vcam1","Gdf15","Slc40a1"))

row.names(markers) <- 1:nrow(markers) 
markers2 <- markers[c(2,3,5,7,8,10,13,15,18,19,23,24,27,28,31,32,34,35,
                      38,40,42,46,50,51,53,55,59,60,61,62,65,66),]

Idents(Myeloid) <- "subcelltype2"
p <- DotPlot(Myeloid,
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
p$data$features.plot <- factor(p$data$features.plot,levels = c("Sell","Mmp9",
                                                               "Icam1","Cxcr2",
                                                               "Hilpda","Vegfa",
                                                               "Cxcl10","Isg20",
                                                               "Xcr1","Btla",
                                                               "Cd209a","Mgl2",
                                                               "Ccr7","Fscn1",
                                                               "Flt3","Tbc1d4",
                                                               "Siglech","Ly6d",
                                                               "Ccr2","Id3",
                                                               "Spp1","Mmp12",
                                                               "Tgfbi","Ocstamp",
                                                               "Ifnb1","Mmp13",
                                                               "Jun","Hspa1a",
                                                               "Fcrls","Vcam1",
                                                               "Gdf15","Slc40a1"))


library(ggh4x)
strip <- strip_themed(background_x = elem_list_rect(fill =c("#A50026","#D73027","#F46D43","#FC4E2A","#FDAE61","#FEE08B","#D9EF8B","#A6D96A","#66BD63","#006837","#9ECAE1","#4292C6","#08306B","#807DBA","#313695","#88419D") ))

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

df <- data.frame(x = 0, y = levels(Myeloid), stringsAsFactors = F )
df$y <- factor(df$y, levels = df$y )

pl <- ggplot(df, aes(x, y, color = factor(y))) +
  geom_point(size = 6, shape = 15,show.legend = F) +
  scale_color_manual(values = rev(c("#A50026","#D73027","#F46D43","#FC4E2A","#FDAE61","#FEE08B","#D9EF8B","#A6D96A","#66BD63","#006837","#9ECAE1","#4292C6","#08306B","#807DBA","#313695","#88419D"
  ))) +
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
p <- plot_grid(pl, p, align = "h", axis="bt", rel_widths = c(1.5, 9))
p

ggsave("dotplot.svg", plot = p, device = "svg",height = 6,width = 13)


#Plot Extended Data Figure 5b
Myeloid <- readRDS("/home/shixi7/zhaochuang/project/pymt/immune/seurat/decontX_immune/Myeloid/Myeloid_final_annotion.rds")
notice <- subset(Myeloid,subcelltype2%in%c("Neu_Sell","Neu_Icam1","Neu_Hilpda_Vegfa","Neu_Cxcl10","Mac_Spp1_Vegfa",
                                           "Mac_Ccr2","Mac_Ifnb1","Mac_Jun","Mac_Fcrls","Mac_Gdf15"))
library(reshape2)
markerdf2=read_tsv("./Marker.tsv")
markerdf2$gene=as.character(markerdf2$gene)
vln.df=as.data.frame(notice[["decontXcounts"]]@data[markerdf2$gene,])
vln.df$gene=rownames(vln.df)
vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")
notice$CB <- colnames(notice)
anno=notice@meta.data[,c("CB","subcelltype2")]
vln.df=inner_join(vln.df,anno,by="CB")
vln.df$gene=factor(vln.df$gene,levels = markerdf2$gene) 

vln.df%>%ggplot(aes(subcelltype2,exp))+geom_violin(aes(fill=gene),scale = "width")+
  facet_grid(vln.df$gene~.,scales = "free_y")+
  scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none"
  )
ggsave("marker_vln.pdf",width = 15,height = 15,units = "cm")



#Plot Extended Data Figure 5c

---
  title: "Survival_Analysis"
author: "Lin Tian"
date: "2024-09-12"
output: html_document
---
  
## Set working directory
  
```{r}
library(survival)
setwd("~/ProTracerDeleter")
load("~/ProTracerDeleter/metabric.exp.RData")
```

## Survival Analysis

```{r}
library(ggplot2)
dim(metabric.a.exp[metabric.a.ann$Pam50Subtype=="Basal", ])
metabric.a.ann$SPP1 <- metabric.a.exp[, "SPP1"]
metabric.a.ann$VEGFA <- metabric.a.exp[, "VEGFA"]
metabric.a.ann$SPP1_g <- ifelse(metabric.a.ann$SPP1 > median(metabric.a.ann$SPP1), "SPP1_high", "SPP1_low")
metabric.a.ann$VEGFA_g <- ifelse(metabric.a.ann$VEGFA > median(metabric.a.ann$VEGFA), "VEGFA_high", "VEGFA_low")
metabric.a.ann$combine_g <- paste(metabric.a.ann$SPP1_g, metabric.a.ann$VEGFA_g, sep="_")


ggplot(metabric.a.ann) + geom_vline(xintercept=median(metabric.a.ann$SPP1)) + geom_hline(yintercept=median(metabric.a.ann$VEGFA)) + geom_point(aes(x=SPP1, y=VEGFA, colour=combine_g)) + scale_color_manual(values=c("#E41A1C", "#386cb0", "#FF7F00", "#4DAF4A")) + theme_bw() + theme(axis.title.x = element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), axis.line=element_line(colour="black"))+coord_fixed(ratio=2)

survdiff(Surv('T'/30, DSS) ~ combine_g, data=metabric.a.ann)
plot(survfit(Surv('T'/30, DSS) ~ combine_g, data=metabric.a.ann), col=c("#E41A1C", "#386cb0", "#FF7F00", "#4DAF4A"), lty=1, xlab = "Months", ylab = "Disease-Free Survival prob.", ylim=c(0, 1), mark.time=TRUE)
legend("topright", c("SPP1_high_VEGFA_high", "SPP1_high_VEGFA_low", "SPP1_low_VEGFA_high", "SPP1_low_VEGFA_low"), col=c("#E41A1C", "#386cb0", "#FF7F00", "#4DAF4A"), lty=1, bty = "n")

```

## Correlation

```{r}
table(metabric.a.ann$Pam50Subtype)

dim(metabric.a.exp[metabric.a.ann$Pam50Subtype=="Basal", ])
metabric.a.ann$SPP1 <- metabric.a.exp[, "SPP1"]
metabric.a.ann$VEGFA <- metabric.a.exp[, "VEGFA"]
metabric.a.ann$SPP1_g <- ifelse(metabric.a.ann$SPP1 > median(metabric.a.ann$SPP1), "SPP1_high", "SPP1_low")
metabric.a.ann$VEGFA_g <- ifelse(metabric.a.ann$VEGFA > median(metabric.a.ann$VEGFA), "VEGFA_high", "VEGFA_low")
metabric.a.ann$combine_g <- paste(metabric.a.ann$SPP1_g, metabric.a.ann$VEGFA_g, sep="_")
survdiff(Surv('T'/30, DSS) ~ combine_g, data=metabric.a.ann)
plot(survfit(Surv('T'/30, DSS) ~ combine_g, data=metabric.a.ann), col=c("#E41A1C", "#386cb0", "#FF7F00", "#4DAF4A"), lty=1, xlab = "Months", ylab = "Disease-Free Survival prob.", ylim=c(0, 1), mark.time=TRUE)
legend("topright", c("SPP1_high_VEGFA_high", "SPP1_high_VEGFA_low", "SPP1_low_VEGFA_high", "SPP1_low_VEGFA_low"), col=c("#E41A1C", "#386cb0", "#FF7F00", "#4DAF4A"), lty=1, bty = "n")

ggplot(metabric.a.ann) + geom_vline(xintercept=median(metabric.a.ann$SPP1)) + geom_hline(yintercept=median(metabric.a.ann$VEGFA)) + geom_point(aes(x=SPP1, y=VEGFA, colour=combine_g)) + scale_color_manual(values=c("#E41A1C", "#386cb0", "#FF7F00", "#4DAF4A")) + theme_bw() + theme(axis.title.x = element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), axis.line=element_line(colour="black"))+coord_fixed(ratio=2)

table(metabric.a.ann$Pam50Subtype)

## Basal
metabric.a.ann.TNBC <- subset(metabric.a.ann, Pam50Subtype=="Basal")
cor.test(metabric.a.ann.TNBC$SPP1, metabric.a.ann.TNBC$VEGFA, method="spearman")
ggplot(metabric.a.ann.TNBC) + geom_point(aes(x=SPP1, y=VEGFA)) + scale_x_continuous(breaks=seq(10, 30, 10), limits=c(10, 30)) + scale_y_continuous(breaks=seq(16, 26, 4), limits=c(16, 26)) + geom_smooth(aes(x=SPP1, y=VEGFA), method="lm") + theme_bw() + theme(axis.title.x = element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), axis.line=element_line(colour="black"))+coord_fixed(ratio=2)

## HER2
metabric.a.ann.HER2 <- subset(metabric.a.ann, Pam50Subtype=="Her2")
cor.test(metabric.a.ann.HER2$SPP1, metabric.a.ann.HER2$VEGFA, method="spearman")
ggplot(metabric.a.ann.HER2) + geom_point(aes(x=SPP1, y=VEGFA)) + scale_x_continuous(breaks=seq(10, 30, 10), limits=c(10, 30)) + scale_y_continuous(breaks=seq(16, 26, 4), limits=c(16, 26)) + geom_smooth(aes(x=SPP1, y=VEGFA), method="lm") + theme_bw() + theme(axis.title.x = element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), axis.line=element_line(colour="black"))+coord_fixed(ratio=2)

## Luminal
metabric.a.ann.Lum <- subset(metabric.a.ann, Pam50Subtype=="LumA" | Pam50Subtype=="LumB")
cor.test(metabric.a.ann.Lum$SPP1, metabric.a.ann.Lum$VEGFA, method="spearman")
ggplot(metabric.a.ann.Lum) + geom_point(aes(x=SPP1, y=VEGFA)) + scale_x_continuous(breaks=seq(10, 30, 10), limits=c(10, 30)) + scale_y_continuous(breaks=seq(16, 26, 4), limits=c(16, 26)) + geom_smooth(aes(x=SPP1, y=VEGFA), method="lm") + theme_bw() + theme(axis.title.x = element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), axis.line=element_line(colour="black"))+coord_fixed(ratio=2)


```




#Plot Extended Data Figure 6a
#Read the RDS
nonimmune <- readRDS("/home/shixi7/zhaochuang/project/pymt/PyMT/normal-annotion/seurat/nonimmune_annotion.rds")
unique(nonimmune$celltype)
Endothelial <- subset(nonimmune,celltype%in%c("Endothelial"))

#Scoring endothelial cells using different signatures
Capillary <- c("Txnip", "Spaca4","Plvap","Socs3")
Capillary = list(Capillary)
Endothelial <- AddModuleScore(Endothelial,
                              features = Capillary,
                              ctrl = 100,
                              name = "Capillary_Features")

Venous <- c("Ccl14", "Ackr1","Postn","Sele","Selp")
Venous = list(Venous)
Endothelial <- AddModuleScore(Endothelial,
                              features = Venous,
                              ctrl = 100,
                              name = "Venous_Features")

Tip <- c("Flt1", "Insr","Pgf","Angpt2")
Tip = list(Tip)
Endothelial <- AddModuleScore(Endothelial,
                              features = Tip,
                              ctrl = 100,
                              name = "Tip_Features")

Immature <- c("Aplnr", "Id1","Eng")
Immature = list(Immature)
Endothelial <- AddModuleScore(Endothelial,
                              features = Immature,
                              ctrl = 100,
                              name = "Immature_Features")

#Compare the Scores between TAM and TAM_DT
library(ggpubr)
for (i in colnames(Endothelial@meta.data)[14:17]) {
  
  use_colors2 <- c("#6BAED6","#D53E4F")
  
  Data = data.frame('intervene'=Endothelial$intervene,
                    'Score'=c(Endothelial@meta.data[,i]))
  library(dplyr)
  Data_summary = Data %>% 
    group_by(intervene) %>% 
    dplyr::summarise(across(.cols = Score,.fns = list(
      'mean' = mean, 'sd' =sd
    ))) %>% 
    dplyr::rename('mean' = Score_mean,
                  'sd' = Score_sd)
  Data = left_join(Data,Data_summary,by="intervene")
  
  compaired <- list(c("TAM", "TAM_DT"))
  
  p1 <- ggplot(Data, aes(x=intervene, y=Score,fill=intervene)) + 
    geom_violin(trim=FALSE,color="white") + #绘制小提琴图
    geom_boxplot(width=0.2,position=position_dodge(0.9))+
    geom_point(aes(x=intervene, y=mean),pch=19,position=position_dodge(0.9),size=1.5)+#绘制均值为点图
    geom_errorbar(aes(ymin = mean-sd, ymax=mean+sd), #误差条表示95%的置信区间
                  width=0.1, #误差条末端短横线的宽度
                  position=position_dodge(0.9), 
                  color="black",
                  alpha = 0.7,
                  size=0.5) +
    scale_fill_manual(values =use_colors2 )+
    geom_signif(comparisons = compaired,
                step_increase = 0.5,
                map_signif_level = F,
                test = "wilcox.test")+#设置填充的颜色
    theme_bw()+ #背景变为白色
    theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times",size=10), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
          panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
          legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                   size=8),
          legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                    size=18),
          panel.grid.major = element_blank(),   #不显示网格线
          panel.grid.minor = element_blank())+  #不显示网格线
    ylab(paste0(i,"_Score"))+xlab("") #设置x轴和y轴的标题
  p1
  
  ggsave(paste0(i,".pdf"), plot = p1, device = "pdf",height = 9,width = 9)
  
}




