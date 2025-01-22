library(tidyverse)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)

library(clusterProfiler)
library(org.Mm.eg.db)

setwd("~/GitHub/MACanalyzeR_Methods/Analysis_RNAseq/GSE111105_PpargKO_Il4_Rosi/")

##################
# LOAD DATA
##################
#LOAD DATA FROM TSV
data <- read.table('countMatrix.txt', header = T)
head(data)
colSums(data[,3:ncol(data)])

#CREATE MATRIX 
countdata <- as.matrix(data[,c(22:25,30:33,26:29,7:10,15:17,11:14)])
countdata <- as.matrix(data[,c(22:25,26:29,7:10,11:14)])

head(countdata)
mode(countdata) <- "integer"
rownames(countdata) <- matrix(unlist(strsplit(data$Geneid, '[.]')),ncol=2,byrow=T)[,1]
rownames(countdata) <- data$gene_name
head(countdata)

#############################
# DIFFERENTIAL GENE EXPRESSION with DESeq2
#############################
# CREATE COLDATA WITH CONDITIONS
condition <- factor(rep(c("WT_DMSO", "WT_Rosi", "WT_IL4", "PPARGko_DMSO", "PPARGko_Rosi", "PPARGko_IL4"), 
                        times =c(4,4,4,4,3,4)))
condition <- factor(rep(c("WT_DMSO", "WT_IL4", "PPARGko_DMSO", "PPARGko_IL4"), 
                        times =c(4,4,4,4)))

# CREATE DESEQ2 OBJ
dds <- DESeqDataSetFromMatrix(
  countData = countdata, 
  colData = as.data.frame(condition), 
  design = ~ condition
)

# SET REFERENCE SAMPLES
dds$condition <- relevel(dds$condition, ref = "WT_DMSO")

# RUN DESEQ2 PLOTS
dds <- DESeq(dds)
resultsNames(dds)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

## NORMALIZED COUNTS
norm <- counts(dds, normalized=TRUE)

############################
## QUALITY CONTROL
############################
# PCA how sample clusters
vsdata <- vst(dds, blind=FALSE)
DESeq2::plotPCA(vsdata)

##################
# RESULTS - DEGs
##################
# wt rosi
res_wt <- results(dds, contrast = c("condition", "WT_Rosi", "WT_DMSO"))
res_wt <- res_wt[which(res_wt$padj < 0.05),]
res_wt <- res_wt[order(res_wt$log2FoldChange, decreasing = T),]

# write.table(res_wt, file="WTrosi_DEG.csv", sep=",", quote=F, col.names=T, row.names = T)

# ko rosi
res_ko <- results(dds, contrast = c("condition", "PPARGko_Rosi", "PPARGko_DMSO"))
res_ko <- res_ko[which(res_ko$padj < 0.05),]
res_ko <- res_ko[order(res_ko$log2FoldChange, decreasing = T),]

# write.table(res_wt, file="PPARGkorosi_DEG.csv", sep=",", quote=F, col.names=T, row.names = T)


## uprgulated
int <- intersect(rownames(res_wt[which(res_wt$log2FoldChange>0),]), rownames(res_ko[which(res_ko$log2FoldChange>0),]))

wt_diff <- setdiff(rownames(res_wt[which(res_wt$log2FoldChange>0),]), rownames(res_ko[which(res_ko$log2FoldChange>0),]))
ko_diff <- setdiff(rownames(res_ko[which(res_ko$log2FoldChange>0),]), rownames(res_wt[which(res_wt$log2FoldChange>0),]))

k_wt <- enrichGO(wt_diff, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "ALL")
barplot(k_wt, showCategory = 10)
View(k_wt@result)

k_ko <- enrichGO(ko_diff, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
barplot(k_ko, showCategory = 10)


# ############# VS LAM
# library(Seurat)
# mac <- readRDS("/home/andrea/Documenti/COLD/elife_dicembre2023.rds")
# lam <- FindMarkers(mac, ident.1 = "LAM")
# lam <- lam[order(lam$avg_log2FC, decreasing = T),]
# head(lam)
# write.table(lam, "../MAC_lam_deg.csv", sep = ",", quote = F, row.names = T, col.names = T)
# 
# int_lam <- intersect(rownames(lam[which(lam$avg_log2FC>0),]), wt_diff)
# k_intlam <- enrichGO(int_lam, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
# barplot(k_intlam, showCategory = 10)
# 
# lam_diff <- setdiff(rownames(lam[which(lam$avg_log2FC>0),]), wt_diff)
# k_difflam <- enrichGO(lam_diff, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
# barplot(k_difflam, showCategory = 10)
# 
# rosi_diff <- setdiff(wt_diff, rownames(lam[which(lam$avg_log2FC>0),]))
# k_diffrosi <- enrichGO(rosi_diff, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
# barplot(k_diffrosi, showCategory = 10)


########################
# GENE SET ENRICHMENT
########################
geneList <- res_wt[,2]
names(geneList) <- rownames(res_wt)
geneList <- sort(geneList, decreasing = T)

ego_all <- gseGO(geneList = geneList, OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL")
ego_bp <- gseGO(geneList = geneList, OrgDb = org.Mm.eg.db, ont = "BP", keyType = "SYMBOL")
ego_mf <- gseGO(geneList = geneList, OrgDb = org.Mm.eg.db, ont = "MF", keyType = "SYMBOL")
ego_cc <- gseGO(geneList = geneList, OrgDb = org.Mm.eg.db, ont = "CC", keyType = "SYMBOL")

dotplot(ego_all, showCategory = 15, font.size = 15, split=".sign", title = "ALL") +
  facet_grid(.~.sign) +
  theme(legend.text = element_text(size = 15), text = element_text(size=25))

dotplot(ego_bp, showCategory = 15, font.size = 15, split=".sign", title = "BP") +
  facet_grid(.~.sign) +
  theme(legend.text = element_text(size = 15), text = element_text(size=25))

dotplot(ego_cc, showCategory = 15, font.size = 15, split=".sign", title = "CC") +
  facet_grid(.~.sign) +
  theme(legend.text = element_text(size = 15), text = element_text(size=25))

dotplot(ego_mf, showCategory = 15, font.size = 15, split=".sign", title = "MF") +
  facet_grid(.~.sign) +
  theme(legend.text = element_text(size = 15), text = element_text(size=25))


res <- res_wt
k_bp <- enrichGO(gene = rownames(res[which(res$log2FoldChange>0),]), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
k_cc <- enrichGO(gene = rownames(res[which(res$log2FoldChange>0),]), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "CC")
k_mf <- enrichGO(gene = rownames(res[which(res$log2FoldChange>0),]), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "MF")

barplot(k_bp, showCategory = 10)
barplot(k_cc,showCategory = 15)
barplot(k_mf,showCategory = 10)


######################
## HEATMAP PLOTTING
######################
#genes
lam_genes <- c("Trem2", "Cd36", "Lpl", "Lipa", "Fabp4", "Plin2", "Pnpla2", "Gdf15")

lam_genes <- c("Ctsa", "Ctsb", "Ctsk", "Ctsl", "Ctss", "Lamp1", "Cd63")

lam_genes <- c("mt-Atp6", "mt-Atp8", "mt-Co1", "mt-Co2", "mt-Co3", "mt-Nd1", "mt-Nd2", "mt-Nd3", "mt-Nd4", "mt-Nd5", "mt-Cytb")


ann_col <- data.frame("Sample" = rep(c("WT_DMSO", "WT_Rosi", "WT_IL4", "PPARGko_DMSO", "PPARGko_Rosi", "PPARGko_IL4"), 
                                     times =c(4,4,4,4,3,4)))
ann_col <- data.frame("Sample" = rep(c("WT_DMSO", "WT_IL4", "PPARGko_DMSO", "PPARGko_IL4"), 
                                     times =c(4,4,4,4)))
rownames(ann_col) <- colnames(countdata)

ann_color <- list("Sample" = c("WT_DMSO" = "#8644A2", "WT_Rosi"="#D862BC", "WT_IL4"="#E59BE9",
                               "PPARGko_DMSO"="#5356FF", "PPARGko_Rosi"="#378CE7", "PPARGko_IL4"="#67C6E3"))

pheatmap(log2(counts(dds, normalized=TRUE)[lam_genes,] + 1), 
         cluster_rows = F, 
         cluster_cols = F, 
         scale = "row",
         show_rownames = T, 
         show_colnames = F,
         annotation_col = ann_col,
         annotation_colors = ann_color, 
         annotation_names_row = F,
         annotation_names_col = F, border_color = F,
         color = colorRampPalette(c("blue", "white", "red"))(100))







