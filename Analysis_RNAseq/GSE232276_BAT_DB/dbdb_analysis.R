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
head(data)

#CREATE MATRIX 
countdata <- as.matrix(data[,3:10])
mode(countdata) <- "integer"

rownames(countdata) <- matrix(unlist(strsplit(data$Geneid, '[.]')),ncol=2,byrow=T)[,1]
rownames(countdata) <- data$gene_name
head(countdata)

#############################
# DIFFERENTIAL GENE EXPRESSION with DESeq2
#############################
# CREATE COLDATA WITH CONDITIONS
condition <- factor(rep(c("WT", "DBDB"), times =c(4,4)))

# CREATE DESEQ2 OBJ
dds <- DESeqDataSetFromMatrix(
  countData = countdata, 
  colData = as.data.frame(condition), 
  design = ~ condition
)

# SET REFERENCE SAMPLES
dds$condition <- relevel(dds$condition, ref = "WT")

# RUN DESEQ2 PLOTS
dds <- DESeq(dds)
resultsNames(dds)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

## NORMALIZED COUNTS
norm <- counts(dds, normalized=TRUE)

############################
## PCA
############################
# PCA how sample clusters
vsdata <- vst(dds, blind=FALSE)
DESeq2::plotPCA(vsdata)

##################
# RESULTS - DEGs
##################
res <- results(dds, contrast = c("condition", "DBDB", "WT"))
res <- res[which(res$padj < 0.05),]
res <- res[order(res$log2FoldChange, decreasing = T),]

write.table(res, file="Foam_DEG.csv", sep=",", quote=F, col.names=T, row.names = T)

########################
# GENE SET ENRICHMENT
########################
geneList <- res[,2]
names(geneList) <- rownames(res)
geneList <- sort(geneList, decreasing = T)

ego_all <- gseGO(geneList = geneList, OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL")
dotplot(ego_all, showCategory = 5, font.size = 15, split=".sign", title = "ALL") +
  facet_grid(.~.sign) +
  theme(legend.text = element_text(size = 15), text = element_text(size=25))

ego_bp <- gseGO(geneList = geneList, OrgDb = org.Mm.eg.db, ont = "BP", keyType = "SYMBOL")
dotplot(ego_bp, showCategory = 5, font.size = 15, split=".sign", title = "BP") +
  facet_grid(.~.sign) +
  theme(legend.text = element_text(size = 15), text = element_text(size=25))

ego_mf <- gseGO(geneList = geneList, OrgDb = org.Mm.eg.db, ont = "MF", keyType = "SYMBOL")
dotplot(ego_cc, showCategory = 5, font.size = 15, split=".sign", title = "CC") +
  facet_grid(.~.sign) +
  theme(legend.text = element_text(size = 15), text = element_text(size=25))

ego_cc <- gseGO(geneList = geneList, OrgDb = org.Mm.eg.db, ont = "CC", keyType = "SYMBOL")
dotplot(ego_mf, showCategory = 5, font.size = 15, split=".sign", title = "MF") +
  facet_grid(.~.sign) +
  theme(legend.text = element_text(size = 15), text = element_text(size=25))


k_bp <- enrichGO(gene = rownames(res[which(res$log2FoldChange>1),]), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
barplot(k_bp, showCategory = 10)

k_cc <- enrichGO(gene = rownames(res[which(res$log2FoldChange>0),]), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "CC")
barplot(k_cc,showCategory = 10)

k_mf <- enrichGO(gene = rownames(res[which(res$log2FoldChange>0),]), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "MF")
barplot(k_mf,showCategory = 10)


## KEGG
k_genes <- mapIds(x=org.Mm.eg.db, keys = rownames(res[which(res$log2FoldChange<0),]), column = "ENTREZID", keytype = "ENSEMBL")
k_genes <- unname(k_genes)

kg <- enrichKEGG(gene = k_genes, organism = "mmu")
kg@result <- separate(kg@result, col = "Description", into = c("Description", "ORG"), sep = " - ")
barplot(kg, showCategory = 5)
View(kg@result)

############# VS LAM of Macrophages single-cell
library(Seurat)

lam <- FindMarkers(mac, ident.1 = "LAM", min.pct = 0.2)
lam <- lam[order(lam$avg_log2FC, decreasing = T),]

int_lam <- intersect(rownames(lam[which(lam$avg_log2FC>1),]), rownames(res[which(res$log2FoldChange>1),]))
k_intlam <- enrichGO(int_lam, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
barplot(k_intlam, showCategory = 5)

lam_diff <- setdiff(rownames(lam[which(lam$avg_log2FC>1),]), rownames(res[which(res$log2FoldChange>1),]))
k_difflam <- enrichGO(lam_diff, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
barplot(k_difflam, showCategory = 7)

foam_diff <- setdiff(rownames(res[which(res$log2FoldChange>1),]), rownames(lam[which(lam$avg_log2FC>1),]))
k_diffrosi <- enrichGO(foam_diff, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")
barplot(k_diffrosi, showCategory = 7)


######################
## HEATMAP PLOTTING
######################
##  brown genes
brown_genes <- c("Lpgat1", "Gria3", "Sgpp1", "Ar", "Alcam", "Hoxc8", "Ndufab1", "Cpt2", 
               "Pdha1", "Shmt1", "Etfa", "Ndufa8", "Ndufs3", "Idh2", "Pgk1", "Cs", "Pald1", 
               "Ndufa9", "Mecr", "Dnajc15", "Ppif", "Fahd1", "Ndufa10", "Idh3g", 
               "Etfb", "Pkm", "Dld", "Uqcr11", "Cycs", "Ndufb8", "Acads", "Fh1", 
               "Sirt3", "Sdhd", "Slc25a19", "Ndufs1", "Cyc1", "Mrpl34", "Ndufv1", 
               "Dlst", "Sdhb", "Etfdh", "Uqcrc1", "Idh3a", "Mrps36", "Hadhb", "Aco2", 
               "Slc25a20", "Impdh1", "Zic1", "Acaa2", "Pdk4", "Cox7a1", "Cidea", "Ucp1")


db <- log2(counts(dds, normalized=TRUE)[brown_genes,] + 1)
db <- counts(dds, normalized=TRUE)[brown_genes,]

ann_db <-  ComplexHeatmap::HeatmapAnnotation(Sample = rep(c("WT", "DB"), times=c(4,4)), show_annotation_name = F, 
                                             col = list(Sample=c("WT"="#416D19", "DB"="#9BCF53")))
ann_genes <-  ComplexHeatmap::rowAnnotation(Genes = rep(c("White Adipose", "Brown Adipose"), times=c(6,49)), show_annotation_name = F,
                                            col = list(Genes=c("White Adipose"="#c7ad7f", "Brown Adipose"="#8B4411")))

ht_db <- ComplexHeatmap::Heatmap(t(scale(t(db))),  name = "norm",
                                   cluster_rows = F, 
                                   cluster_columns = F, 
                                   show_column_names = F, 
                                   top_annotation = ann_db, right_annotation = ann_genes)

ComplexHeatmap::draw(ht_db)



##################################################àà
##  PathwayScore Analysis
##################################################àà
pathway_path <- '../pathways/'
path <- list.files(path=pathway_path, recursive = F, full.names = F)
path <- matrix(unlist(str_split(path, "[.]")), ncol = 2, byrow = T)[,1]

V1 <- matrix(NA, nrow=length(colnames(norm)), ncol=length(path), dimnames=list(colnames(norm), path))
norm <- norm[which(rowSums(norm)>0),]

for (p in path) {
  genes <- as.data.frame(read_tsv(paste0(pathway_path, p, ".txt")))[,2]
  genes <- intersect(genes, rownames(norm))
  
  sample_mean <- apply(mtx_pathway, 1, function(x)by(x, cond_tot, mean))
  #ratio_expr <- t(sample_mean) / colMeans(sample_mean)
  gene_weight <- apply(mtx_pathway, 1, var)
  mean_exp_pathway <- apply(mtx_pathway, 2, function(x) weighted.mean(x, gene_weight/sum(gene_weight)))
  V1[,p] <-  mean_exp_pathway
  }

##  correlation
cors_db <- cor(V1, method = "pearson")

col_fun = circlize::colorRamp2(c(-1, 0, 1), c("darkgreen", "white", "purple4"))
anncor_db <-  ComplexHeatmap::HeatmapAnnotation(Sample = rep(c("Overfeeding"), times=c(4)), show_annotation_name = F,
                                                 col = list(Sample=c("Overfeeding"="#416D19")))
cor_db <- ComplexHeatmap::Heatmap(cors_db, name = "correlation", col = col_fun, rect_gp = gpar(type = "none"), 
                                     cell_fun = function(j, i, x, y, width, height, fill) {
                                       grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = NA, fill = NA))
                                       grid.circle(x = x, y = y, r = abs(cors_db[i, j])/2 * min(unit.c(width, height)), 
                                                   gp = gpar(fill = col_fun(cors_db[i, j]), col = NA))
                                     }, 
                                     cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left",
                                     show_row_names = T, show_column_names = T, column_names_side = "top", bottom_annotation = anncor_db)

ComplexHeatmap::draw(cor_db)

##  linear model
fit_db <- lm(Macrophage_Activation ~ Adaptive_Thermogenesis, data = V1)
drop1(fit_db, test = "Chisq")
performance::check_model(fit_db)
visreg::visreg(fit_db, xvar = "Adaptive_Thermogenesis")


