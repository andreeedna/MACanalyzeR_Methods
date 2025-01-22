library(tidyverse)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)

library(clusterProfiler)
library(org.Mm.eg.db)

setwd("~/GitHub/MACanalyzeR_MaterialsMethods/Analysis_RNAseq/GSE_PpargKOrosi/")

##################
# LOAD DATA
##################
#LOAD DATA FROM TSV
data <- read.table('countMatrix.txt', header = T)
head(data)
colSums(data[,3:ncol(data)])

#CREATE MATRIX 
countdata <- as.matrix(data[,c(22:25,30:33,7:10,15:17)])
head(countdata)
mode(countdata) <- "integer"
rownames(countdata) <- matrix(unlist(strsplit(data$Geneid, '[.]')),ncol=2,byrow=T)[,1]
rownames(countdata) <- data$gene_name
head(countdata)

#############################
# DIFFERENTIAL GENE EXPRESSION with DESeq2
#############################
# CREATE COLDATA WITH CONDITIONS
condition <- factor(rep(c("WT_DMSO", "WT_Rosi", "PPARGko_DMSO", "PPARGko_Rosi"), times =c(4,4,4,3)))

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
# write.table(normalized_counts, file="countMatrix_normalized.txt", sep="\t", quote=F, col.names=NA)

##################
# RESULTS - DEGs
##################
# wt rosi
res_wt <- results(dds, contrast = c("condition", "WT_Rosi", "WT_DMSO"))
res_wt <- res_wt[which(res_wt$padj < 0.05),]
res_wt <- res_wt[order(res_wt$log2FoldChange, decreasing = T),]

# write.table(res_wt, file="WTrosi_DEG.csv", sep=",", quote=F, col.names=T, row.names = T)

# ko rosi
res_ko <- results(dds, contrast = c("condition", "PPARGko_Rosi", "WT_DMSO"))
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
lam_genes <- c("Trem2", "Cd36", "Lpl", "Lipa", "Fabp4", "Plin2", "Pnpla2")
lam_genes <- c("Ctsa", "Ctsb", "Ctsk", "Ctsl", "Ctss", "Lamp1", "Cd63")

lam_genes <- c("mt-Atp6", "mt-Atp8", "mt-Co1", "mt-Co2", "mt-Co3", "mt-Nd1", "mt-Nd2", "mt-Nd3", "mt-Nd4", "mt-Nd5", "mt-Cytb")


ann_col <- data.frame("Sample" = rep(c("WT_DMSO", "WT_Rosi", "PPARGko_DMSO", "PPARGko_Rosi"), times =c(4,4,4,3)))
rownames(ann_col) <- colnames(countdata)

ann_color <- list("Sample" = c("WT_DMSO" = "#210B2C", "WT_Rosi"="#55286F", "PPARGko_DMSO"="#BC96E6", "PPARGko_Rosi"="#D8B4E2"))

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


#############

genes <- MACanalyzeR::mouse_go$LIPID_METABOLIC_PROCESS
genes <- intersect(genes, rownames(countdata))

ann_col <- data.frame("Sample" = rep(c("WT_DMSO", "WT_Rosi", "PPARGko_DMSO", "PPARGko_Rosi"), times =c(4,4,4,3)))
rownames(ann_col) <- colnames(countdata)

ann_color <- list("Sample" = c("WT_DMSO" = "#210B2C", "WT_Rosi"="#55286F", "PPARGko_DMSO"="#BC96E6", "PPARGko_Rosi"="#D8B4E2"))


pmat <- log2(counts(dds, normalized=TRUE)[genes,] + 1)
#pmat <- counts(dds, normalized=TRUE)[genes,]

pmat <- pmat[is.finite(rowSums(pmat)),]
pmat <- pmat[which(rowSums(pmat)!=0),]

pheatmap(pmat, 
         cluster_rows = T, 
         cluster_cols = F, 
         scale = "row",
         show_rownames = T, 
         show_colnames = F,
         annotation_col = ann_col,
         annotation_colors = ann_color, 
         annotation_names_row = F,
         annotation_names_col = F, 
         color = colorRampPalette(c("blue", "white", "red"))(100))




# path
lyso_gen <- matrix(unlist(strsplit(read.table("lysosomal_genes.txt", header = F, sep = "\t")[,2], ";")), ncol = 2, byrow = T)[,1]
lyso_gen <- intersect(lyso_gen, rownames(countdata))

oxphos_gen <- matrix(unlist(strsplit(read.table("oxphos_genes.txt", header = F, sep = "\t")[,2], ";")), ncol = 2, byrow = T)[,1]
oxphos_gen <- intersect(oxphos_gen, rownames(countdata))

lipid_gen <-  matrix(unlist(strsplit(read.table("lipiddegr_genes.txt", header = F, sep = "\t")[,2], ";")), ncol = 2, byrow = T)[,1]
lipid_gen <- intersect(lipid_gen, rownames(countdata))

lyso_gen <- setdiff(setdiff(lyso_gen, lipid_gen), oxphos_gen) 
lipid_gen <- setdiff(setdiff(lipid_gen, lyso_gen), oxphos_gen) 
oxphos_gen <- setdiff(setdiff(oxphos_gen, lipid_gen), lyso_gen) 

gen <- c(lipid_gen, lyso_gen, oxphos_gen)

mexpr <- log2(counts(dds, normalized=TRUE)[gen,] + 1)
mexpr <- mexpr[which(rowSums(mexpr)!=0),]

ann_row <- data.frame("Sample" = rep(c("WT_DMSO", "WT_Rosi", "PPARGko_DMSO", "PPARGko_Rosi"), times =c(4,4,4,3)))
rownames(ann_row) <- colnames(countdata)

t <- c(length(intersect(rownames(mexpr), lipid_gen)), 
       length(intersect(rownames(mexpr), lyso_gen)), 
       length(intersect(rownames(mexpr), oxphos_gen)))

ann_col <- data.frame("Pathway" = rep(c("Lipid", "Lysosome", "OXPHOS"), times = t))
rownames(ann_col) <- rownames(mexpr)

ann_color <- list("Sample" = c("WT_DMSO" = "#210B2C", "WT_Rosi"="#55286F", "PPARGko_DMSO"="#BC96E6", "PPARGko_Rosi"="#D8B4E2"),
                  "Pathway" = c("Lipid"="gold", "Lysosome"="darkgreen", "OXPHOS"="orange"))

pheatmap::pheatmap(mexpr, scale = "row", 
                   annotation_col = ann_row, annotation_names_col = F,
                   annotation_row = ann_col, annotation_names_row = F,
                   annotation_colors = ann_color, annotation_legend = T,
                   cluster_rows = F, cluster_cols = F, show_colnames = F, 
                   show_rownames = F, border_color = "none")




