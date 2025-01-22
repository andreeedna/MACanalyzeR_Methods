library(tidyverse)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(grid)

library(clusterProfiler)
library(org.Mm.eg.db)

setwd("~/GitHub/MACanalyzeR_MaterialsMethods/Analysis_RNAseq/GSE112740_BAT_HFD/")

##################
# LOAD DATA
##################
#LOAD DATA FROM TSV
data <- read.table('countMatrix.txt', header = T)
head(data)
colSums(data[,3:ncol(data)])

#CREATE MATRIX 
countdata <- as.matrix(data[,c(6:8,3:5)])
head(countdata)
mode(countdata) <- "integer"
rownames(countdata) <- matrix(unlist(strsplit(data$Geneid, '[.]')),ncol=2,byrow=T)[,1]
rownames(countdata) <- data$gene_name
head(countdata)

countdata <- countdata[- grep("^Gm[1234567890]", rownames(countdata)),]
countdata <- countdata[- grep("Rik", rownames(countdata)),]
countdata <- countdata[- grep("^ENSMUSG", rownames(countdata)),]
head(countdata)

#############################
# DIFFERENTIAL GENE EXPRESSION with DESeq2
#############################
# CREATE COLDATA WITH CONDITIONS
condition <- factor(rep(c("LFD", "HFD"), times =c(3,3)))

# CREATE DESEQ2 OBJ
dds <- DESeqDataSetFromMatrix(
  countData = countdata, 
  colData = as.data.frame(condition), 
  design = ~ condition
)

# SET REFERENCE SAMPLES
dds$condition <- relevel(dds$condition, ref = "LFD")

# RUN DESEQ2 PLOTS
dds <- DESeq(dds)
resultsNames(dds)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

## NORMALIZED COUNTS
norm <- counts(dds, normalized=TRUE)
# write.table(normalized_counts, file="countMatrix_normalized.txt", sep="\t", quote=F, col.names=NA)

############################
## QUALITY CONTROL
############################
# PCA how sample clusters
vsdata <- vst(dds, blind=FALSE)
DESeq2::plotPCA(vsdata)

##################
# RESULTS - DEGs
##################
res <- results(dds, contrast = c("condition", "HFD", "LFD"))
res <- res[which(res$padj < 0.05),] # 0.049178
res <- res[order(res$log2FoldChange, decreasing = T),]

openxlsx::write.xlsx(as.data.frame(res), "hfd_deg.xlsx", rowNames=T)

# ############# VS LAM
# library(Seurat)
# mac <- readRDS("../db_EFSD/DB_HFD_secondaversione/Macrophages.rds")
# DimPlot(mac)
# lam <- FindMarkers(mac, ident.1 = "LAM", min.pct = 0.2)
# lam <- lam[order(lam$avg_log2FC, decreasing = T),]
# head(lam)
# 
# int_lam <- intersect(rownames(lam[which(lam$avg_log2FC>1),]), rownames(res[which(res$log2FoldChange>1),]))
# k_intlam <- enrichGO(int_lam, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
# k_intlam@result[c(1,6,8,9),] <- NA
# barplot(k_intlam, showCategory = 5)
# View(k_intlam@result)
# 
# lam_diff <- setdiff(rownames(lam[which(lam$avg_log2FC>1),]), rownames(res[which(res$log2FoldChange>1),]))
# k_difflam <- enrichGO(lam_diff, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
# k_difflam@result[c(4,6:10),] <- NA
# barplot(k_difflam, showCategory = 7)
# 
# 
# foam_diff <- setdiff(rownames(res[which(res$log2FoldChange>1),]), rownames(lam[which(lam$avg_log2FC>1),]))
# k_diffrosi <- enrichGO(foam_diff, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
# barplot(k_diffrosi, showCategory = 7)


########################
# GENE SET ENRICHMENT
########################
geneList <- res[,2]
names(geneList) <- rownames(res)
geneList <- sort(geneList, decreasing = T)

ego_all <- gseGO(geneList = geneList, OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL")
ego_bp <- gseGO(geneList = geneList, OrgDb = org.Mm.eg.db, ont = "BP", keyType = "SYMBOL")
ego_mf <- gseGO(geneList = geneList, OrgDb = org.Mm.eg.db, ont = "MF", keyType = "SYMBOL")
ego_cc <- gseGO(geneList = geneList, OrgDb = org.Mm.eg.db, ont = "CC", keyType = "SYMBOL")

dotplot(ego_all, showCategory = 5, font.size = 15, split=".sign", title = "ALL") +
  facet_grid(.~.sign) +
  theme(legend.text = element_text(size = 15), text = element_text(size=25))

dotplot(ego_bp, showCategory = 5, font.size = 15, split=".sign", title = "BP") +
  facet_grid(.~.sign) +
  theme(legend.text = element_text(size = 15), text = element_text(size=25))

dotplot(ego_cc, showCategory = 5, font.size = 15, split=".sign", title = "CC") +
  facet_grid(.~.sign) +
  theme(legend.text = element_text(size = 15), text = element_text(size=25))

dotplot(ego_mf, showCategory = 5, font.size = 15, split=".sign", title = "MF") +
  facet_grid(.~.sign) +
  theme(legend.text = element_text(size = 15), text = element_text(size=25))


k_bp <- enrichGO(gene = rownames(res[which(res$log2FoldChange>1),]), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
k_cc <- enrichGO(gene = rownames(res[which(res$log2FoldChange>0),]), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "CC")
k_mf <- enrichGO(gene = rownames(res[which(res$log2FoldChange>0),]), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "MF")

barplot(k_bp, showCategory = 10)
barplot(k_cc,showCategory = 10)
barplot(k_mf,showCategory = 10)



## KEGG
k_genes <- mapIds(x=org.Mm.eg.db, keys = rownames(res[which(res$log2FoldChange>0),]), column = "ENTREZID", keytype = "ENSEMBL")
k_genes <- unname(k_genes)

kg <- enrichKEGG(gene = k_genes, organism = "mmu")
kg@result <- separate(kg@result, col = "Description", into = c("Description", "ORG"), sep = " - ")
#kg@result[c(2,3,5:11, 13:16),] <- NA
barplot(kg, showCategory = 20)
View(kg@result)


######################
## HEATMAP PLOTTING
######################
#genes
lam_genes <- c("Trem2", "Cd36", "Cd9", "Lpl", "Lipa", "Fabp4", "Plin2", "Pnpla2", "Spp1")
lam_genes <- c("Ctsa", "Ctsb", "Ctsk", "Ctsl", "Ctss", "Lamp1", "Cd63")
lam_genes <- c("mt-Atp6", "mt-Atp8", "mt-Co1", "mt-Co2", "mt-Co3", "mt-Nd1", "mt-Nd2", "mt-Nd3", "mt-Nd4", "mt-Nd5", "mt-Cytb")

ann_col <- data.frame("Sample" = rep(c("LFD", "HFD"), times =c(3,3)))
rownames(ann_col) <- colnames(countdata)

#ann_color <- list("Sample" = c("WT_DMSO" = "#210B2C", "WT_Rosi"="#55286F", "PPARGko_DMSO"="#BC96E6", "PPARGko_Rosi"="#D8B4E2"))

pheatmap(log2(counts(dds, normalized=TRUE)[lam_genes,] + 1), 
         cluster_rows = T, 
         cluster_cols = F, 
         scale = "row",
         show_rownames = T, 
         show_colnames = F,
         annotation_col = ann_col,
         #annotation_colors = ann_color, 
         annotation_names_row = F,
         annotation_names_col = F, 
         color = colorRampPalette(c("blue", "white", "red"))(100))


######################
## HEATMAP PLOTTING
######################
lam_genes <- c("Lpgat1", "Gria3", "Sgpp1", "Ar", "Alcam", "Hoxc8", "Ndufab1", "Cpt2", 
               "Pdha1", "Shmt1", "Etfa", "Ndufa8", "Ndufs3", "Idh2", "Pgk1", "Cs", "Pald1", 
               "Ndufa9", "Mecr", "Dnajc15", "Ppif", "Fahd1", "Ndufa10", "Idh3g", 
               "Etfb", "Pkm", "Dld", "Uqcr11", "Cycs", "Ndufb8", "Acads", "Fh1", 
               "Sirt3", "Sdhd", "Slc25a19", "Ndufs1", "Cyc1", "Mrpl34", "Ndufv1", 
               "Dlst", "Sdhb", "Etfdh", "Uqcrc1", "Idh3a", "Mrps36", "Hadhb", "Aco2", 
               "Slc25a20", "Impdh1", "Zic1", "Acaa2", "Pdk4", "Cox7a1", "Cidea", "Ucp1")


hfd <- log2(counts(dds, normalized=TRUE)[lam_genes,c(6,2,3,4,5,1)] + 1)
hfd <- counts(dds, normalized=TRUE)[lam_genes,]

ann_hfd <-  ComplexHeatmap::HeatmapAnnotation(Sample = rep(c("LFD", "HFD"), times=c(3,3)), show_annotation_name = F,
                                              col = list(Sample=c("LFD"="#E25E3E", "HFD"="#FF9B50")))

ht_hfd <- ComplexHeatmap::Heatmap(t(scale(t(hfd))),  name = "norm",
                                  cluster_rows = F, 
                                  cluster_columns = F, 
                                  show_column_names = T, 
                                  top_annotation = ann_hfd)

ComplexHeatmap::draw(ht_hfd)

##################################################àà
pathway_path <- '../pathways/prova/'
path <- list.files(path=pathway_path, recursive = F, full.names = F)
path <- matrix(unlist(str_split(path, "[.]")), ncol = 2, byrow = T)[,1]

V1 <- matrix(NA, nrow=length(colnames(norm)), ncol=length(path), dimnames=list(colnames(norm), path))
norm <- norm[which(rowSums(norm)>0),]

for (p in path) {
  genes <- as.data.frame(read_tsv(paste0(pathway_path, p, ".txt")))[,2]
  genes <- intersect(genes, rownames(norm))
  V1[,p] <- colMeans(norm[genes,])
}

cors_db <- cor(V1, method = "pearson")

col_fun = circlize::colorRamp2(c(-1, 0, 1), c("darkgreen", "white", "purple4"))
# anncor_db <-  ComplexHeatmap::HeatmapAnnotation(Sample = rep(c("Overfeeding"), times=c(4)), show_annotation_name = F,
#                                                  col = list(Sample=c("Overfeeding"="#416D19")))
cor_db <- ComplexHeatmap::Heatmap(cors_db, name = "correlation", col = col_fun, rect_gp = gpar(type = "none"), 
                                     cell_fun = function(j, i, x, y, width, height, fill) {
                                       grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = NA, fill = NA))
                                       grid.circle(x = x, y = y, r = abs(cors_db[i, j])/2 * min(unit.c(width, height)), 
                                                   gp = gpar(fill = col_fun(cors_db[i, j]), col = NA))
                                     }, 
                                     cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left",
                                     show_row_names = T, show_column_names = T, column_names_side = "top", #bottom_annotation = anncor_db
                                  )

ComplexHeatmap::draw(cor_db)

###

























