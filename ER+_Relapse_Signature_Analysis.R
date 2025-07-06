---
title: "ER_relapse"
author: "Lin Tian"
date: "2025-04-26"
output: html_document
editor_options: 
chunk_output_type: console
---
  
## Load Data
  
  ```{r}
library(genefilter)
library(ggplot2)
setwd("~/Tamoxifen")
load("~/Tamoxifen/GSE12665.RData")
all.equal(colnames(GSE12665.exp), rownames(GSE12665.ann))
```

## PCA plot

```{r}
GSE12665.ann_sub <- subset(GSE12665.ann, ER==1)

rowVariance_combine <- rowVars(GSE12665.exp[, rownames(GSE12665.ann_sub)])
names(rowVariance_combine) <- rownames(GSE12665.exp_sub)
selected_Rows_combine <- names(sort(rowVariance_combine, decreasing = T)[1:500])
cx_combine <- sweep(t(GSE12665.exp[selected_Rows_combine, rownames(GSE12665.ann_sub)]), 2, colMeans(t(GSE12665.exp[selected_Rows_combine, rownames(GSE12665.ann_sub)])))
sv_combine <- svd(cx_combine)
sv_combine$d^2/sum(sv_combine$d^2) # variance explained
GSE12665.ann_sub$PC1 <- sv_combine$u[,1]
GSE12665.ann_sub$PC2 <- sv_combine$u[,2]

ggplot(GSE12665.ann_sub) + geom_point(aes(x=PC1, y=PC2, colour=Treatment), size=6) + scale_color_manual(values=c("#AC0000","#57A5E5", "#999999"))
```

## Signature

```{r}
GSE12665.ann_sub <- subset(GSE12665.ann, ER==1)
GSE12665.ann_sub$ESR1 <- GSE12665.exp["ESR1", rownames(GSE12665.ann_sub)]
t.test(GSE12665.ann_sub$ESR1~ GSE12665.ann_sub$Treatment)
Esr1_plot <- ggplot(data=GSE12665.ann_sub, aes(x=Treatment, y=ESR1)) + geom_boxplot(fill=NA, outliers = FALSE) + geom_point(aes(colour=Treatment), size=6, shape=16, alpha=0.5) + scale_color_manual(values=c("#57A5E5", "#AC0000")) + scale_y_continuous(limits=c(6, 14), breaks=seq(6, 14, by=4)) +  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), axis.line.x=element_line(colour="black"), axis.line.y=element_line(colour="black")) + coord_fixed(ratio=0.5)
ggsave("~/Esr1_plot.pdf", Esr1_plot, width=6.28, height=6.28, units="in", useDingbats=FALSE)


CD8_effector_sig <- c("CD8A", "CD8B", "GZMA", "GZMB", "GZMK", "PRF1", "IFNG", "GNLY", "CCL5", "EOMES", "TBX21", "FASLG", "ZAP70")
CD8_effector_sig <- CD8_effector_sig[CD8_effector_sig %in% rownames(GSE12665.exp)]
GSE12665.ann_sub$CD8_effector <- apply(GSE12665.exp[CD8_effector_sig, rownames(GSE12665.ann_sub)], 2, sum)
t.test(GSE12665.ann_sub$CD8_effector ~ GSE12665.ann_sub$Treatment)
CD8_effector_plot <- ggplot(data=GSE12665.ann_sub, aes(x=Treatment, y=CD8_effector)) + geom_boxplot(fill=NA, outliers = FALSE) + geom_point(aes(colour=Treatment), size=6, shape=16, alpha=0.5) + scale_color_manual(values=c("#57A5E5", "#AC0000")) + scale_y_continuous(limits=c(60, 100), breaks=seq(60, 100, by=10)) +  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), axis.line.x=element_line(colour="black"), axis.line.y=element_line(colour="black")) + coord_fixed(ratio=0.1)
ggsave("~/CD8_effector_plot.pdf", CD8_effector_plot, width=6.28, height=6.28, units="in", useDingbats=FALSE)

CD8_exhausted_sig <- c("CD8A", "CD8B", "PDCD1", "LAYN", "HAVCR2", "LAG3", "CD244", "CTLA4", "LILRB1", "TIGIT", "TOX", "VSIR", "BTLA", "ENTPD1", "LAIR1")
CD8_exhausted_sig <- CD8_exhausted_sig[CD8_exhausted_sig %in% rownames(GSE12665.exp)]
GSE12665.ann_sub$CD8_exhausted <- apply(GSE12665.exp[CD8_exhausted_sig, rownames(GSE12665.ann_sub)], 2, sum)
t.test(GSE12665.ann_sub$CD8_exhausted ~ GSE12665.ann_sub$Treatment)
CD8_exhausted_plot <- ggplot(data=GSE12665.ann_sub, aes(x=Treatment, y=CD8_exhausted)) + geom_boxplot(fill=NA, outliers = FALSE) + geom_point(aes(colour=Treatment), size=6, shape=16, alpha=0.5) + scale_color_manual(values=c("#57A5E5", "#AC0000")) + scale_y_continuous(limits=c(60, 100), breaks=seq(60, 100, by=10)) +  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), axis.line.x=element_line(colour="black"), axis.line.y=element_line(colour="black")) + coord_fixed(ratio=0.1)
ggsave("~/CD8_exhausted_plot.pdf", CD8_exhausted_plot, width=6.28, height=6.28, units="in", useDingbats=FALSE)

NK_sig <- c("NCR1", "NCAM1", "FCGR3A", "NKG7", "GNLY", "KLRD1", "KLRF1", "KLRC1", "KLRC2")
NK_sig <- NK_sig[NK_sig %in% rownames(GSE12665.exp)]
GSE12665.ann_sub$NK <- apply(GSE12665.exp[NK_sig, rownames(GSE12665.ann_sub)], 2, sum)
t.test(GSE12665.ann_sub$NK ~ GSE12665.ann_sub$Treatment)
NK_plot <- ggplot(data=GSE12665.ann_sub, aes(x=Treatment, y=NK)) + geom_point() + geom_boxplot() + geom_boxplot(fill=NA, outliers = FALSE) + geom_point(aes(colour=Treatment), size=6, shape=16, alpha=0.5) + scale_color_manual(values=c("#57A5E5", "#AC0000")) + scale_y_continuous(limits=c(40, 70), breaks=seq(40, 70, by=10)) +  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), axis.line.x=element_line(colour="black"), axis.line.y=element_line(colour="black")) + coord_fixed(ratio=0.13)
ggsave("~/NK_plot.pdf", NK_plot, width=6.28, height=6.28, units="in", useDingbats=FALSE)

gdT_sig <- c("TRDC", "TRDV1", "TRDV2", "TRDV3", "TRGC1", "TRGC2", "TRGV9", "CD3D", "CD3E", "CD3G")
gdT_sig <- gdT_sig[gdT_sig %in% rownames(GSE12665.exp)]
GSE12665.ann_sub$gdT <- apply(GSE12665.exp[gdT_sig, rownames(GSE12665.ann_sub)], 2, sum)
t.test(GSE12665.ann_sub$CD8_exhausted ~ GSE12665.ann_sub$Treatment)
gdT_plot <- ggplot(data=GSE12665.ann_sub, aes(x=Treatment, y=gdT)) + geom_boxplot() + geom_boxplot(fill=NA, outliers = FALSE) + geom_point(aes(colour=Treatment), size=6, shape=16, alpha=0.5) + scale_color_manual(values=c("#57A5E5", "#AC0000")) + scale_y_continuous(limits=c(8, 20), breaks=seq(8, 20, by=4)) +  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), axis.line.x=element_line(colour="black"), axis.line.y=element_line(colour="black")) + coord_fixed(ratio=0.33)
ggsave("~/gdT_plot.pdf", gdT_plot, width=6.28, height=6.28, units="in", useDingbats=FALSE)
```