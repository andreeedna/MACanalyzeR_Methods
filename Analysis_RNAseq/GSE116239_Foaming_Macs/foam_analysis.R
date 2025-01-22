library(DESeq2)

library(tidyverse)
library(openxlsx)
library(ggplot2)
library(ggrepel)

library(clusterProfiler)
library(org.Mm.eg.db)

##################
# LOAD DATA
##################
#LOAD DATA FROM TSV
data <- read.table('countMatrix.txt', header = T)
colSums(data[,3:ncol(data)])
head(data)

#CREATE MATRIX 
countdata <- as.matrix(data[,c(6:8,3:5)])
head(countdata)
mode(countdata) <- "integer"
# rownames(countdata) <- matrix(unlist(strsplit(data$Geneid, '[.]')),ncol=2,byrow=T)[,1]
rownames(countdata) <- data$gene_name
head(countdata)

#############################
# DIFFERENTIAL GENE EXPRESSION with DESeq2
#############################
# CREATE COLDATA WITH CONDITIONS
condition <- factor(rep(c("NonFoam", "Foam"), times =c(3,3)))

# CREATE DESEQ2 OBJ
dds <- DESeqDataSetFromMatrix(
  countData = countdata, 
  colData = as.data.frame(condition), 
  design = ~ condition
)

# SET REFERENCE SAMPLES
dds$condition <- relevel(dds$condition, ref = "NonFoam")

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
res <- results(dds, contrast = c("condition", "Foam", "NonFoam"))
res <- res[which(res$padj < 0.05),]
res <- res[order(res$log2FoldChange, decreasing = T),]

write.table(res, file="Foam_DEG.csv", sep=",", quote=F, col.names=T, row.names = T)

############# VS LAM of Macrophages single-cell
library(Seurat)
library(MACanalyzeR)

mac_efsd <- CreateMacObj(mac, "Sample", "Cells")
mac_efsd <- FoamSpotteR(mac_efsd)
mac$Foam <- mac_efsd@MetaData$Foam

DimPlot(mac, group.by = "Foam")

lam <- FindMarkers(mac, ident.1 = "fMAC+", group.by = "Foam", min.pct = 0.4, subset.ident = "LAM")
lam <- lam[order(lam$avg_log2FC, decreasing = T),]
head(lam)

int_lam <- intersect(rownames(lam[which(lam$avg_log2FC>0.25),]), rownames(res[which(res$log2FoldChange>1),]))
k_intlam <- enrichGO(int_lam, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
barplot(k_intlam, showCategory = 7, font.size = 17) + 
  theme(legend.text = element_text(size = 17), legend.title = element_text(size=17))

lam_diff <- setdiff(rownames(lam[which(lam$avg_log2FC>0.25),]), rownames(res[which(res$log2FoldChange>1),]))
k_difflam <- enrichGO(lam_diff, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
barplot(k_difflam, showCategory = 7, font.size = 17) + 
  theme(legend.text = element_text(size = 17), legend.title = element_text(size=17))


foam_diff <- setdiff(rownames(res[which(res$log2FoldChange>0.25),]), rownames(lam[which(lam$avg_log2FC>1),]))
k_diffrosi <- enrichGO(foam_diff, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
barplot(k_diffrosi, showCategory = 7, font.size = 17) + 
  theme(legend.text = element_text(size = 17), legend.title = element_text(size=17))

