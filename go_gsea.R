# ------------------------------------------------------------------
# Single-cell RNA-seq GO and GSEA Enrichment Analysis
# ------------------------------------------------------------------

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(tidyverse)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(clusterProfiler) # GO, KEGG, GSEA analysis
  library(enrichplot)      # Visualization for enrichment results
  library(GSEABase)        # GSEA analysis
  library(GSVA)            # Gene Set Variation Analysis
  library(DOSE)            # Disease Ontology Semantic and Enrichment analysis
  library(topGO)           # GO enrichment
  library(msigdbr)         # MSigDB gene sets
  library(pheatmap)        # Heatmaps
})

# Load Seurat object for epithelial cells
Epithelial <- readRDS("/home/shixi7/zhaochuang/project/pymt/PyMT/normal-annotion/seurat/Epithelial/Epithelial_annotion.rds")
unique(Epithelial$subcelltype1)

# Define cell groups: group 1 = "Cancer_cells_c10_Ly6a", group 2 = all others
cells1 <- subset(Epithelial@meta.data, subcelltype1 == "Cancer_cells_c10_Ly6a") %>% rownames()
cells2 <- subset(Epithelial@meta.data, subcelltype1 != "Cancer_cells_c10_Ly6a") %>% rownames()

# Find differentially expressed genes (DEGs) between two groups
deg <- FindMarkers(Epithelial, ident.1 = cells1, ident.2 = cells2)
deg <- data.frame(gene = rownames(deg), deg)
head(deg)

# Filter DEGs by adjusted p-value and log fold change thresholds
deg1 <- deg
k1 <- (deg1$p_val_adj < 0.05) & (deg1$avg_log2FC < -0.25)  # Downregulated genes
k2 <- (deg1$p_val_adj < 0.05) & (deg1$avg_log2FC > 0.25)   # Upregulated genes
table(k1)
table(k2)

# Assign gene regulation status
deg1$change <- ifelse(k1, "down", ifelse(k2, "up", "stable"))
head(deg1)

# Convert gene symbols to Entrez IDs (mouse)
s2e <- bitr(deg1$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
deg1 <- inner_join(deg1, s2e, by = c("gene" = "SYMBOL"))
head(deg1)

# Prepare gene lists for GO and KEGG analysis
gene_up       <- deg1[deg1$change == "up", "gene"]
gene_down     <- deg1[deg1$change == "down", "gene"]
gene_diff     <- c(gene_up, gene_down)

gene_all       <- deg1[, "ENTREZID"]
gene_up_KEGG   <- deg1[deg1$change == "up", "ENTREZID"]
gene_down_KEGG <- deg1[deg1$change == "down", "ENTREZID"]
gene_diff_KEGG <- c(gene_up_KEGG, gene_down_KEGG)

# GO enrichment analysis (Cellular Component, upregulated genes)
ego_CC <- enrichGO(gene         = gene_up,
                   keyType      = "SYMBOL",
                   OrgDb        = org.Mm.eg.db,
                   ont          = "CC",
                   pAdjustMethod= "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)
head(data.frame(ego_CC))

# GO enrichment analysis (Biological Process, upregulated genes)
ego_BP <- enrichGO(gene         = gene_up,
                   keyType      = "SYMBOL",
                   OrgDb        = org.Mm.eg.db,
                   ont          = "BP",
                   pAdjustMethod= "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)
head(data.frame(ego_BP))

# GO enrichment analysis (Molecular Function, upregulated genes)
ego_MF <- enrichGO(gene         = gene_up,
                   keyType      = "SYMBOL",
                   OrgDb        = org.Mm.eg.db,
                   ont          = "MF",
                   pAdjustMethod= "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)

# Combined GO enrichment for all ontologies (up- and downregulated)
go     <- enrichGO(gene = gene_up, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "ALL")
godown <- enrichGO(gene = gene_down, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "ALL")

# Dotplot visualization split by ontology
p <- dotplot(go, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scales = "free")
print(p)

p <- dotplot(godown, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scales = "free")
print(p)

# Barplot for biological processes (upregulated)
p <- barplot(ego_BP, showCategory = 20, colorBy = "p.adjust")
print(p)
ggsave(p, file = "Ly6a_GO_BP.svg", width = 10, height = 10)

# Select specific GO terms to highlight
select_terms <- c("extracellular matrix organization",
                  "extracellular structure organization",
                  "cell-matrix adhesion",
                  "mesenchyme development",
                  "epithelial cell migration",
                  "mesenchymal cell differentiation",
                  "connective tissue development")

# Filter enrichment results for selected terms and order by adjusted p-value
ego_select <- ego_BP@result[ego_BP@result$Description %in% select_terms, ]
ego_select <- ego_select[order(ego_select$p.adjust), ]

# Replace GO results with filtered subset for plotting
ego_select_obj <- ego_BP
ego_select_obj@result <- ego_select

# Plot selected GO terms as barplot
p <- barplot(ego_select_obj, showCategory = nrow(ego_select), title = "Selected GO Terms")
print(p)


# ----------------------------------GSEA--------------------------------

# Load mouse Hallmark gene sets (category "H") from msigdbr
genesets <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  as.data.frame()

# Check available species and gene set categories
msigdbr_collections()

# Prepare named ranked gene list for GSEA (log2 fold changes sorted descending)
geneLists <- deg1$avg_log2FC
names(geneLists) <- deg1$gene
geneLists <- sort(geneLists, decreasing = TRUE)

# Perform GSEA using pre-defined gene sets, with p-value cutoff 0.1
egmt <- GSEA(geneLists, TERM2GENE = genesets, verbose = FALSE, pvalueCutoff = 0.1)

# Filter significantly enriched pathways (adjusted p-value < 0.25)
significant_pathways <- egmt@result %>% filter(p.adjust < 0.25)
head(significant_pathways)

# Convert GSEA results to dataframe for inspection
gsea_df <- data.frame(egmt)
head(gsea_df)

# Plot classic GSEA enrichment plot for pathway with ID = 5
p <- gseaplot2(egmt, geneSetID = 5, pvalue_table = TRUE)
print(p)
ggsave(p, file = "Ly6a_GSEA.svg", width = 11, height = 10)

# Filter positively enriched pathways (NES > 0) with significance
positive_pathways <- egmt@result %>% filter(NES > 0, p.adjust < 0.25)
head(positive_pathways)

# Plot GSEA enrichment for the second positively enriched pathway
p <- gseaplot2(egmt, geneSetID = positive_pathways$ID[2], pvalue_table = TRUE)
print(p)
ggsave(p, file = "Ly6a_GSEA.svg", width = 8, height = 8)

