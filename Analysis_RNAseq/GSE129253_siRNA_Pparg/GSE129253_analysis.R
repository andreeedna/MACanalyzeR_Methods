library(tidyverse)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)

library(clusterProfiler)
library(org.Mm.eg.db)

setwd("~/GitHub/MACanalyzeR_Methods/Analysis_RNAseq/GSE129253_siRNA_Pparg/")

##################
# LOAD DATA
##################
#LOAD DATA FROM TSV
data <- read.table('GSE129253_Pparg_matrix.txt', header = T)
head(data)
colSums(data[,3:ncol(data)])

#CREATE MATRIX 
countdata <- as.matrix(data[,c(3:11)])
countdata <- as.matrix(data[,c(3:5, 9:11)])
head(countdata)
mode(countdata) <- "integer"
rownames(countdata) <- data$gene_name
head(countdata)

# countdata <- countdata[- grep("^Gm[1234567890]", rownames(countdata)),]
# countdata <- countdata[- grep("Rik", rownames(countdata)),]

#############################
# DIFFERENTIAL GENE EXPRESSION with DESeq2
#############################
# CREATE COLDATA WITH CONDITIONS
condition <- factor(rep(c("CTRL", "CTRL_neg", "Pparg_siPool"), times =c(3,3,3)))
condition <- factor(rep(c("CTRL", "Pparg_siPool"), times =c(3,3)))

# CREATE DESEQ2 OBJ
dds <- DESeqDataSetFromMatrix(
  countData = countdata, 
  colData = as.data.frame(condition), 
  design = ~ condition
)

# SET REFERENCE SAMPLES
dds$condition <- relevel(dds$condition, ref = "CTRL")

# RUN DESEQ2 PLOTS
dds <- DESeq(dds)
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
res <- results(dds, contrast = c("condition", "CTRL", "Pparg_siPool"))
res <- res[which(res$padj < 0.05),]
res <- res[order(res$log2FoldChange, decreasing = T),]

# write.table(res, file="ctrl_vs_pparg.csv", sep=",", quote=F, col.names=T, row.names = T)

## cluster profiler
k_u <- enrichGO(gene = rownames(res[which(res$log2FoldChange>0),]), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "ALL")
barplot(k_u, showCategory=10)

######################
## HEATMAP PLOTTING
######################
#genes
lam_genes <- c("Trem2", "Cd36", "Cd9", "Lpl", "Lipa", "Fabp4", "Plin2", "Pnpla2", "Gdf15")
lam_genes <- c("Ctsa", "Ctsb", "Ctsk", "Ctsl", "Ctss", "Lamp1", "Cd63")
lam_genes <- c("mt-Atp6", "mt-Atp8", "mt-Co1", "mt-Co2", "mt-Co3", "mt-Nd1", "mt-Nd2", "mt-Nd3", "mt-Nd4", "mt-Nd5", "mt-Cytb")

lam_genes <- intersect(lam_genes, rownames(countdata))

# ann_col <- data.frame("Sample" = rep(c("ctrl", "ctrl_neg", "siPparg"), times =c(3,3,3)))
ann_col <- data.frame("Sample" = rep(c("ctrl", "siRNA Pparg"), times =c(3,3)))
rownames(ann_col) <- colnames(countdata)

ann_color <- list("Sample" = c("ctrl"="#850F8D", "siRNA Pparg"="#C738BD"))

pheatmap(log2(counts(dds, normalized=TRUE)[lam_genes,] + 1), 
         cluster_rows = F, 
         cluster_cols = F, 
         scale = "row",
         show_rownames = T, 
         show_colnames = F,
         annotation_col = ann_col,
         annotation_colors = ann_color, 
         annotation_names_row = F, fontsize = 15,
         annotation_names_col = F, border_color = F,
         color = colorRampPalette(c("blue", "white", "red"))(100))



########################àà
# path
lyso_gen <- matrix(unlist(strsplit(read.table("../GSE_PpargKOrosi/lysosomal_genes.txt", header = F, sep = "\t")[,2], ";")), ncol = 2, byrow = T)[,1]
lyso_gen <- intersect(lyso_gen, rownames(countdata))

oxphos_gen <- matrix(unlist(strsplit(read.table("../GSE_PpargKOrosi/oxphos_genes.txt", header = F, sep = "\t")[,2], ";")), ncol = 2, byrow = T)[,1]
oxphos_gen <- intersect(oxphos_gen, rownames(countdata))

lipid_gen <-  matrix(unlist(strsplit(read.table("../GSE_PpargKOrosi/lipiddegr_genes.txt", header = F, sep = "\t")[,2], ";")), ncol = 2, byrow = T)[,1]
lipid_gen <- intersect(lipid_gen, rownames(countdata))

lyso_gen <- setdiff(setdiff(lyso_gen, lipid_gen), oxphos_gen) 
lipid_gen <- setdiff(setdiff(lipid_gen, lyso_gen), oxphos_gen) 
oxphos_gen <- setdiff(setdiff(oxphos_gen, lipid_gen), lyso_gen) 

gen <- c(lipid_gen, lyso_gen, oxphos_gen)

mexpr <- log2(counts(dds, normalized=TRUE)[gen,c(1:3,7:9)] + 1)
mexpr <- mexpr[which(rowSums(mexpr)!=0),]

ann_row <- data.frame("Sample" = rep(c("ctrl", "ctrl_neg", "siPparg"), times =c(3,3,3)))
rownames(ann_row) <- colnames(countdata)

t <- c(length(intersect(rownames(mexpr), lipid_gen)), 
       length(intersect(rownames(mexpr), lyso_gen)), 
       length(intersect(rownames(mexpr), oxphos_gen)))

ann_col <- data.frame("Pathway" = rep(c("Lipid", "Lysosome", "OXPHOS"), times = t))
rownames(ann_col) <- rownames(mexpr)

ann_color <- list(# "Sample" = c("WT_DMSO" = "#210B2C", "WT_Rosi"="#55286F", "PPARGko_DMSO"="#BC96E6", "PPARGko_Rosi"="#D8B4E2"),
                  "Pathway" = c("Lipid"="gold", "Lysosome"="darkgreen", "OXPHOS"="orange"))

pheatmap::pheatmap(mexpr, scale = "row", 
                   annotation_col = ann_row, 
                   annotation_names_col = F,
                   annotation_row = ann_col, annotation_names_row = F,
                   annotation_colors = ann_color, annotation_legend = T,
                   cluster_rows = F, cluster_cols = F, show_colnames = F, 
                   show_rownames = F, border_color = "none")




