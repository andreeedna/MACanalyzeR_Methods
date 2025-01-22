#scRNAseq BAT tissue integrated analysis
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

set.seed(500)

tot_color_sample <-c("wt"="#E69F00", "hfd"="#56B4E9", "db/db"="#009E73")
mac_color <- c("Monocytes"="#3468C0", "PVM"="#86A7FC", "LAM"="#FF9843", "P_LAM"="#FFDD95")

##########################
# LOAD DATA
##########################
#get data location
dirs <- list.dirs(path='rawdata', recursive = F, full.names = F)

for(x in dirs){
  name <- gsub('','',x)
  cts <- ReadMtx(mtx = paste0("rawdata/", x, '/filtered_feature_bc_matrix/matrix.mtx.gz'),
                 features = paste0("rawdata/", x,'/filtered_feature_bc_matrix/features.tsv.gz'),
                 cells = paste0("rawdata/", x, '/filtered_feature_bc_matrix/barcodes.tsv.gz'))
  assign(x, CreateSeuratObject(counts = cts, project = x))
  rm(cts, name, x)
}

rm(dirs)

###############################
## QC and FILTERING
###############################
# calculate % MT genes
wt[['percent.mt']] <- PercentageFeatureSet(wt, pattern = '^mt-')
hfd[['percent.mt']] <- PercentageFeatureSet(hfd, pattern = '^mt-')
dbdb16w[['percent.mt']] <- PercentageFeatureSet(dbdb16w, pattern = '^mt-')

# filter wt
wt_filt <- subset(wt, subset= nFeature_RNA > 250 & nCount_RNA > 500 & nCount_RNA < 30000 & percent.mt < 10)
VlnPlot(wt_filt, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
wt_filt

# filter hfd
hfd_filt <- subset(hfd, subset= nFeature_RNA > 250 & nCount_RNA > 500 & nCount_RNA < 30000 & percent.mt < 10)
VlnPlot(hfd_filt, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
hfd_filt

# filter dbdb16w
dbdb16w_filt <- subset(dbdb16w, subset= nFeature_RNA > 250 & nCount_RNA > 500 & nCount_RNA < 30000 & percent.mt < 10)
VlnPlot(dbdb16w_filt, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
dbdb16w_filt


# sample the same number of cells - wt16w
n.cell.sample <- length(rownames(hfd_filt@meta.data))
n.cell.sample

sampled.cells <- sample(rownames(wt_filt@meta.data), size = n.cell.sample, replace = F)
wt_filt <- subset(wt_filt, cells = sampled.cells)
wt_filt

sampled.cells <- sample(rownames(dbdb16w_filt@meta.data), size = n.cell.sample, replace = F)
dbdb16w_filt <- subset(dbdb16w_filt, cells = sampled.cells)
dbdb16w_filt

rm(n.cell.sample, sampled.cells)
rm(wt, hfd, dbdb16w)


#############################
# MERGE SEURAT OBJECTS
#############################
adipo <- merge(x= wt_filt , y= c(dbdb16w_filt, hfd_filt),
                       add.cell.id = c('wt', 'db/db', "hfd"))

adipo$sample <- rownames(adipo@meta.data)
adipo@meta.data <- separate(adipo@meta.data, col='sample',
                                    into=c('Sample','Barcode'),
                                    sep='_')

rm(wt_filt, dbdb16w_filt, hfd_filt)

#############################
# INTEGRATION STEPS
#############################
adipo <- NormalizeData(adipo)
adipo <- FindVariableFeatures(adipo)
adipo <- ScaleData(adipo)
adipo <- RunPCA(adipo)
adipo <- IntegrateLayers(
  object = adipo, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = T)

ElbowPlot(adipo, reduction = "pca", 50)

adipo <- FindNeighbors(adipo, reduction = "integrated.cca", dims = 1:20)
adipo <- FindClusters(adipo, resolution = 0.5)
adipo <- RunUMAP(adipo, reduction = "integrated.cca", dims = 1:20)

adipo <- JoinLayers(adipo)

DimPlot(adipo,
        group.by = c("Sample", "seurat_clusters"),
        combine = T, label.size = 2
)

adipo$Sample <- factor(adipo$Sample, levels = c("wt", "hfd", "db/db"))
DimPlot(adipo, split.by = "Sample")

adipo <- subset(adipo, idents = 9, invert=T)
adipo <- FindClusters(adipo, resolution = 0.5)
DimPlot(adipo, label = T)
DimPlot(adipo, label = T, split.by = "Sample")


#################################
## SEARCH AND FIND ALL MARKERS
#################################
totmarker <- c("Adgre1", "Itgam", "Lyz2", #mono/MACs
               "Dpp4", "Sema3c", "Scara5", #ASC
               "Lum", "Podn", "Col5a3", #FAP
               'Trbc2', 'Cd3d', 'Cd3g', #Tcells
              'Igkc', 'Ighm', 'Iglc2', #B cells
              'S100a8', 'S100a9', 'Retnlg', #Neutrophils
              "Mki67", 'Birc5', 'Top2a', #MSC
               'Cldn5', 'Cdh5', 'Pecam1' #endothelials
)

pmat <- apply(adipo@assays$RNA$data[totmarker,], 1, function(x)by(x,adipo$seurat_clusters,mean))
pheatmap::pheatmap(t(pmat), scale = "row", cluster_rows = F, cluster_cols = T, color = colorRampPalette(c("blue", "white", "red"))(100))

#set cluster names
new.cluster.id <- c('FAP', 'Macrophages', 'Endothelial Cells', "ASC", "Macrophages", "ASC", "Macrophages",
                    "Macrophages", "ASC", "T Cells", "Endothelial Cells", "FAP")
adipo$Cells <- factor(adipo$seurat_clusters, levels = 0:11, labels = new.cluster.id)
Idents(adipo) <- "Cells"
rm(new.cluster.id)
DimPlot(adipo)

#save the object
saveRDS(adipo, 'annotato.rds')

################
# MACROPHAGE
################
mac <- subset(adipo, idents = c("Macrophages"))
mac <- NormalizeData(mac)
mac <- FindVariableFeatures(mac)
mac <- ScaleData(mac)
mac <- RunPCA(mac)
mac <- IntegrateLayers(
  object = mac, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = T)

ElbowPlot(mac, reduction = "pca", 50)

mac <- FindNeighbors(mac, reduction = "integrated.cca", dims = 1:10)
mac <- FindClusters(mac, resolution = 0.2)
mac <- RunUMAP(mac, reduction = "integrated.cca", dims = 1:10)

mac <- JoinLayers(mac)

DimPlot(mac,
        group.by = c("Sample", "seurat_clusters"),
        combine = T, label.size = 2
)

# # mac markers
# VlnPlot(mac, c("Adgre1", "Itgam", "Lyz2", "Trem2", "Cd36"), pt.size = 0, ncol = 2)
# 
# # mono markers
# VlnPlot(mac, c('Ccr2', 'Cd44', 'Cx3cr1'),  pt.size = 0)
# FeaturePlot(mac, c('Ccr2', 'Cd44', 'Cx3cr1'),  pt.size = 0)

#macrophage subpopulation
sub <- c("Mrc1", "Lyve1", "Cd163" #PVM
         ,"Cd36", "Trem2", "Cd9" #LAM 
         ,"Pola1", "Kif11", "Kif15", #P-LAM
         'Ccr2', 'Cd44', 'Cx3cr1' #mono
)

pmat <- apply(mac@assays$RNA$data[sub,], 1, function(x)by(x,mac$seurat_clusters,mean))
pheatmap::pheatmap(t(pmat), scale = "row", cluster_rows = F, cluster_cols = T, color = colorRampPalette(c("blue", "white", "red"))(100), 
                   border_color = F)


clu <- c("LAM", "PVM", "Monocytes", "P_LAM")
mac$Cells <- factor(mac$seurat_clusters, levels = 0:3, labels = clu)
Idents(mac) <- "Cells"
DimPlot(mac)

#save the object
saveRDS(mac, 'Macrophages.rds')

#########################################################
##  VELOCITY ANALYSIS
#########################################################
# SCVELO prep
mac$barcode <- colnames(mac)
mac$UMAP_1 <- mac@reductions$umap@cell.embeddings[,1]
mac$UMAP_2 <- mac@reductions$umap@cell.embeddings[,2]
write.csv(mac@meta.data, file='MAC_velocity/metadata_tot.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(mac, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0('velocyto/counts_tot.mtx'))

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(mac@reductions$pca@cell.embeddings, file='velocyto/pca_tot.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='velocyto/gene_names_tot.csv',
  quote=F,row.names=F,col.names=F
)

########################################################
## MACanalyzeR analysis
########################################################
library(MACanalyzeR)
mac_efsd <- CreateMacObj(mac, "Sample", "Cells")
MacPlot(mac_efsd, col = mac_color)


# FoamSpotteR
mac_efsd <- FoamSpotteR(mac_efsd)
FoamPlot(mac_efsd)
FoamPlot(mac_efsd, split.by = "Sample", ncol = 3)


# MacPolarizeR
mac_efsd <- MacPolarizeR(mac_efsd)
MacRadar(mac_efsd, plot.by = "Sample")
MacPlot(mac_efsd, plot.by = "Mac", split.by = "Sample")

dat <- as.data.frame(mac_efsd@MacPolarizeR$SingleCell)
dat[,"Sample"] <- mac_efsd@MetaData[["Sample"]]
dat[,"Cluster"] <- mac_efsd@MetaData[["Cluster"]]

ggplot(dat, aes(x=dat[,1], y=dat[,2])) +
  labs(x="M1", y="M2") +
  geom_vline(xintercept = 1, colour = 'grey30') +
  geom_hline(yintercept = 1, colour = 'grey30') +
  theme(text = element_text(size = 17),
        panel.background = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.line = element_line(),
        strip.background = element_rect(color = "black")) +
  geom_point(aes(color = Sample, fill = Sample), size=3) +
  scale_x_continuous(limits = c(-0.1, 2.5)) +
  scale_y_continuous(limits = c(0, 2.5)) +
  scale_color_manual(values = tot_color_sample) +
  facet_wrap(~Sample)
  

ggplot(dat, aes(x=dat[,1], y=dat[,2])) +
  labs(x="M1", y="M2") +
  geom_vline(xintercept = 1, colour = 'grey30') +
  geom_hline(yintercept = 1, colour = 'grey30') +
  theme(text = element_text(size = 17),
        panel.background = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.line = element_line(),
        strip.background = element_rect(color = "black")) +
  geom_point(aes(color = Cluster, fill = Cluster), size=3) +
  scale_x_continuous(limits = c(-0.1, 2.5)) +
  scale_y_continuous(limits = c(0, 2.5)) +
  scale_color_manual(values = mac_color) +
  facet_grid(rows = vars(Sample), cols = vars(Cluster))


# PathAnalyzeR
#sample
k <- kegg_selected[c(1,2,3,4,5,17)]
k$MITOCHONDRION <- go$MITOCHONDRION
mac_efsd <- PathAnalyzeR(mac_efsd, meta = "Sample", pathway = k)
PathHeat(mac_efsd, plot.by = "Sample", pval = 1)

PathDisplay(mac_efsd)
PathCart(mac_efsd, c(1,5), fill = T, plot.by = "Sample", col = tot_color_sample, a = 0.07)
PathCart(mac_efsd, c(4,2), fill = T, plot.by = "Sample", col = tot_color_sample, a = 0.07)
PathCart(mac_efsd, c(6,7), fill = T, plot.by = "Sample", col = tot_color_sample, a = 0.07)

#cluster
mac_efsd <- PathAnalyzeR(mac_efsd, meta = "Cluster", pathway = k)
PathHeat(mac_efsd, plot.by = "Cluster", pval = 1)

PathDisplay(mac_efsd, meta = "Cluster")
PathCart(mac_efsd, c(1,5), fill = T, plot.by = "Cluster", col = mac_color, a = 0.2) + facet_wrap(~Cluster)
PathCart(mac_efsd, c(4,6), fill = T, plot.by = "Cluster", col = mac_color, a = 0.2) + facet_wrap(~Cluster)

PathPlot(mac_efsd, 1, split.by = "Sample", ncol = 3, max.cutoff = 2.5)
PathPlot(mac_efsd, 4, split.by = "Sample", ncol = 3, max.cutoff = 2.5)
PathPlot(mac_efsd, 5, split.by = "Sample", ncol = 3)
PathPlot(mac_efsd, 6, split.by = "Sample", ncol = 3, max.cutoff = 2.5)

