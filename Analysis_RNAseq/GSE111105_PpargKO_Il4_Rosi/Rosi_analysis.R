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
condition <- factor(rep(c("WT_DMSO", "WT_Rosi", "PPARGko_DMSO", "PPARGko_Rosi"), 
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
#res_wt <- res_wt[which(res_wt$padj < 0.05),]
res_wt <- res_wt[order(res_wt$log2FoldChange, decreasing = T),]

# ko rosi
res_ko <- results(dds, contrast = c("condition", "PPARGko_Rosi", "PPARGko_DMSO"))
#res_ko <- res_ko[which(res_ko$padj < 0.05),]
res_ko <- res_ko[order(res_ko$log2FoldChange, decreasing = T),]



# fc_mat <- cbind(res_wt["Gdf15","log2FoldChange"], res_ko["Gdf15","log2FoldChange"])
# rownames(fc_mat) <- "Gdf15"
# colnames(fc_mat) <- c("WT_Rosi", "PPARGko_Rosi")
# 
# pval_mat <- cbind(res_wt["Gdf15","padj"], res_ko["Gdf15","padj"])
# star_mat <- matrix(cut(pval_mat, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", " ")), ncol = 6, byrow = F)
# 
# gg_mat <- reshape2::melt(fc_mat)
# gg_mat$pval <- reshape2::melt(pval_mat)[,"value"]
# gg_mat$stars <- cut(gg_mat$pval, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
# 
# ggplot(aes(x=Var2, y=Var1, fill=value), data=gg_mat) + 
#   geom_tile() + scale_fill_gradient2(high="#D7191C", mid="white", low="#2C7BB6") + 
#   #   geom_text(aes(label=stars, color=value), size=8) + scale_colour_gradient(low="grey30", high="white", guide="none") +
#   geom_text(aes(label=stars), color="black", size=5) + 
#   labs(y=NULL, x=NULL, fill="log2FC") +
#   theme_minimal() + theme(axis.text.x=element_text(angle = -45, hjust = 0, size=12), axis.text.y=element_text(size=12))
# 
# ## ComplexHeatmap
# # colori
# cold_color <- c("WT_Rosi"="#D862BC", "PPARGko_Rosi"="#378CE7")
# 
# col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#2C7BB6", "white", "#D7191C"))
# 
# # annotazione
# sample_ann <- ComplexHeatmap::HeatmapAnnotation(Sample = c("WT_Rosi", "PPARGko_Rosi"),show_annotation_name = F, 
#                                                 col = list(Sample=cold_color), show_legend = F)
# 
# # legend
# lgd_list = list(
#   ## log2FC
#   ComplexHeatmap::Legend(col_fun = col_fun, title = "log2FC", break_dist = 1,
#                          legend_height = unit(4, "cm"), grid_width = unit(0.5, "cm"),
#                          title_gp = grid::gpar(fontsize=12), labels_gp = grid::gpar(fontsize=13)),
#   ## Sample
#   ComplexHeatmap::Legend(labels = c("WT_Rosi", "PPARGko_Rosi"), title = "Sample", legend_height = unit(0.5, "cm"), grid_width = unit(0.5, "cm"), 
#                          legend_gp = grid::gpar(fill = cold_color), border = F, title_gp = grid::gpar(fontsize=12), 
#                          labels_gp = grid::gpar(fontsize=13))
# )
# 
# fc_heatmap <- ComplexHeatmap::Heatmap(fc_mat, name = "log2FC", col = col_fun,
#                                       cluster_rows = T, cluster_columns = FALSE, 
#                                       row_names_side = "left", column_names_side = "bottom",
#                                       show_row_names = T, show_column_names = F, show_heatmap_legend = F,
#                                       top_annotation = sample_ann,
#                                       cell_fun = function(j, i, x, y, width, height, fill) {
#                                         grid::grid.text(star_mat[i, j], x, y, gp = grid::gpar(fontsize = 10))}
# )
# 
# 
# ComplexHeatmap::draw(fc_heatmap, annotation_legend_list = lgd_list)
# 


######################
## HEATMAP PLOTTING
######################
#genes
lam_genes <- c("Trem2", "Cd36", "Lpl", "Lipa", "Fabp4", "Plin2", "Pnpla2", "Gdf15")
lam_genes <- c("Gdf15")

ann_col <- data.frame("Sample" = rep(c("WT_DMSO", "WT_Rosi", "PPARGko_DMSO", "PPARGko_Rosi"), times =c(4,4,4,4)))
rownames(ann_col) <- colnames(countdata)

ann_color <- list("Sample" = c("WT_DMSO" = "#8644A2", "WT_Rosi"="#D862BC", "PPARGko_DMSO"="#5356FF", "PPARGko_Rosi"="#378CE7"))

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

################################################################################






