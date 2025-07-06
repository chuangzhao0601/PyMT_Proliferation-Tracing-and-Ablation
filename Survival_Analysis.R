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
