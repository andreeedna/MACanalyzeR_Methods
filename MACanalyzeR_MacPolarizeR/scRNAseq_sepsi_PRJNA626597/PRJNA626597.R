## scRNAseq analysis of SVF extracted from SEPTIC adipose tissue after one day and one month
## analysis performed with Seurat 4
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
set.seed(500)

setwd("MACanalyzeR_MacPolarizeR/scRNAseq_sepsi_PRJNA626597/")

######################
# LOAD DATA
######################
#load data
for(x in list.dirs(path='.', recursive = F, full.names = F)){
  name <- gsub('','',x)
  cts <- ReadMtx(mtx = paste0(x, '/matrix.mtx.gz'),
                 features = paste0(x, '/features.tsv.gz'),
                 cells = paste0(x, '/barcodes.tsv.gz'))
  assign(x, CreateSeuratObject(counts = cts, project = x))
  rm(cts, name, x)
}

###############################
# QC and FILTERING
###############################
#calculate % MT genes
CTRL[['percent.mt']] <- PercentageFeatureSet(CTRL, pattern = '^[Mm][Tt]-')
SEPSIS_ONE_DAY[['percent.mt']] <- PercentageFeatureSet(SEPSIS_ONE_DAY, pattern = '^[Mm][Tt]-')
SEPSIS_ONE_MONTH[['percent.mt']] <- PercentageFeatureSet(SEPSIS_ONE_MONTH, pattern = '^[Mm][Tt]-')

# CTRL
VlnPlot(CTRL, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
#filter CTRL
CTRL_filt <- subset(CTRL, subset= nFeature_RNA > 500 & nCount_RNA > 500 & percent.mt < 5)
VlnPlot(CTRL_filt, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
CTRL_filt


# SEPSIS_ONE_DAY
VlnPlot(SEPSIS_ONE_DAY, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
#filter SEPSIS_ONE_MONTH
SEPSIS_ONE_DAY_filt <- subset(SEPSIS_ONE_DAY, subset= nFeature_RNA > 500 & nCount_RNA > 500 & percent.mt < 5)
VlnPlot(SEPSIS_ONE_DAY_filt, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
SEPSIS_ONE_DAY_filt


# SEPSIS_ONE_MONTH
VlnPlot(SEPSIS_ONE_MONTH, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
#filter SEPSIS_ONE_MONTH
SEPSIS_ONE_MONTH_filt <- subset(SEPSIS_ONE_MONTH, subset= nFeature_RNA > 500 & nCount_RNA > 500 & percent.mt < 5)
VlnPlot(SEPSIS_ONE_MONTH_filt, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
SEPSIS_ONE_MONTH_filt


#sample the same number of cells - campione minore WT16w
n.cell.sample <- 2300

sampled.cells <- sample(rownames(CTRL_filt@meta.data), size = n.cell.sample, replace = F)
CTRL_filt <- subset(CTRL_filt, cells = sampled.cells)
CTRL_filt

sampled.cells <- sample(rownames(SEPSIS_ONE_DAY_filt@meta.data), size = n.cell.sample, replace = F)
SEPSIS_ONE_DAY_filt <- subset(SEPSIS_ONE_DAY_filt, cells = sampled.cells)
SEPSIS_ONE_DAY_filt

sampled.cells <- sample(rownames(SEPSIS_ONE_MONTH_filt@meta.data), size = n.cell.sample, replace = F)
SEPSIS_ONE_MONTH_filt <- subset(SEPSIS_ONE_MONTH_filt, cells = sampled.cells)
SEPSIS_ONE_MONTH_filt

rm(n.cell.sample, sampled.cells)
rm(CTRL, SEPSIS_ONE_DAY, SEPSIS_ONE_MONTH)

#############################
# MERGE SEURAT OBJECTS
#############################
sep <- merge(x= CTRL_filt , y=c(SEPSIS_ONE_DAY_filt, SEPSIS_ONE_MONTH_filt),
                       add.cell.id = c('CTRL', 'SEPDAY', 'SEPMONTH'))

sep$sample <- rownames(sep@meta.data)
sep@meta.data <- separate(sep@meta.data, col='sample',
                                    into=c('Sample','Barcode'),
                                    sep='_')

sep$Sample <- factor(x = sep$Sample, levels = c('CTRL', 'SEPDAY', 'SEPMONTH'))

rm(CTRL_filt, SEPSIS_ONE_DAY_filt, SEPSIS_ONE_MONTH_filt)

############################
# INTEGRATION STEP
############################
list <- SplitObject(sep, split.by = 'Sample')
rm(sep)
for(i in 1:length(list)){
  list[[i]] <- NormalizeData(object = list[[i]])
  list[[i]] <- FindVariableFeatures(object = list[[i]], nfeatures = 2000)
  
  rm(i)
}

features <- SelectIntegrationFeatures(object.list = list, nfeatures = 2000)
anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features)
sepsi <- IntegrateData(anchorset = anchors)
sepsi <- ScaleData(object = sepsi)

rm(list ,anchors, features)

DefaultAssay(sepsi) <- "integrated"
sepsi <- ScaleData(object = sepsi)
sepsi <- RunPCA(object = sepsi, npcs = 30)

ElbowPlot(sepsi, reduction = 'pca', 30)

sepsi <- RunUMAP(object = sepsi, reduction = 'pca', dims= 1:20)
sepsi <- FindNeighbors(sepsi, reduction ='pca' , dims = 1:20)
sepsi <- FindClusters(sepsi, resolution = 0.4)

DefaultAssay(sepsi) <- "RNA"

#visualize
DimPlot(sepsi, label = T)
DimPlot(sepsi, split.by = "Sample")

DimPlot(sepsi, label = T)


##  marker annotation
FeaturePlot(sepsi, c("Ptprc", "Pecam1", "Pdgfra"),
            order = T, min.cutoff = "q5", max.cutoff = "q90")

FeaturePlot(sepsi, c("Adgre1", "Mrc1", "Cd163"), order = T, 
            min.cutoff = "q5", max.cutoff = "q90")
FeaturePlot(sepsi, c("Ly6c2", "Ccr2"), order = T, 
            min.cutoff = "q5", max.cutoff = "q90")
FeaturePlot(sepsi, c("Ms4a1", "Cd79a"), order = T, 
            min.cutoff = "q5", max.cutoff = "q90")
FeaturePlot(sepsi, c("Trbc1", "Cd3d", "Nkg7"), order = T, 
            min.cutoff = "q5", max.cutoff = "q90")

FeaturePlot(sepsi, c("Ifi205", "Naaa"), order = T, 
            min.cutoff = "q5", max.cutoff = "q90")
FeaturePlot(sepsi, c("Napsa", "Cd209a"), order = T, 
            min.cutoff = "q5", max.cutoff = "q90")
FeaturePlot(sepsi, c("Siglech", "Ccr9"), order = T, 
            min.cutoff = "q5", max.cutoff = "q90")

FeaturePlot(sepsi, c("S100a8", "S100a9"), order = T, 
            min.cutoff = "q5", max.cutoff = "q90")

FeaturePlot(sepsi, c("Pecam1", "Car4"), order = T, 
            min.cutoff = "q5", max.cutoff = "q90")

new.id <- c("MonoMacs", "Neutrophils", "Fibroblasts", "MonoMacs", "T Cells", "NK", "Neutrophils", "cDC", "T Cells", "cDC", 
            "Neutrophils", "B Cells", "Endothelial Cells", "pDC", "Neurons", "Proliferating")
sepsi$Cells <- factor(sepsi$seurat_clusters, levels = 0:15, labels = new.id)
Idents(sepsi) <- "Cells"

DimPlot(sepsi)

###############
## MONOMAC
###############
mac <- subset(sepsi, idents = c(0,2)) ## clusters expressing Adgre1

list <- SplitObject(mac, split.by = 'Sample')
for(i in 1:length(list)){
  list[[i]] <- NormalizeData(object = list[[i]])
  list[[i]] <- FindVariableFeatures(object = list[[i]], nfeatures = 2000)
  rm(i)
}

features <- SelectIntegrationFeatures(object.list = list, nfeatures = 2000)
anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features)
mac_sepsi <- IntegrateData(anchorset = anchors)
mac_sepsi <- ScaleData(object = mac_sepsi)

rm(list ,anchors, features)

DefaultAssay(mac_sepsi) <- "integrated"
mac_sepsi <- ScaleData(object = mac_sepsi)
mac_sepsi <- RunPCA(object = mac_sepsi, npcs = 30)

ElbowPlot(mac_sepsi, reduction = 'pca', 30)

mac_sepsi <- RunUMAP(object = mac_sepsi, reduction = 'pca', dims= 1:10)
mac_sepsi <- FindNeighbors(mac_sepsi, reduction ='pca' , dims = 1:10)
mac_sepsi <- FindClusters(mac_sepsi, resolution = 0.1)

DimPlot(mac_sepsi)
DimPlot(mac_sepsi, split.by = "Sample")

FeaturePlot(mac_sepsi, "Trem2", order = T, split.by = "Sample", min.cutoff = "q5", max.cutoff = "q90")

##################
##  MACanalyzeR
##################
library(MACanalyzeR)
mac_sep <- CreateMacObj(mac_sepsi, "Sample")
MacPlot(mac_sep, split.by = "Sample")
mac_sep <- MacPolarizeR(mac_sep)
MacPlot(mac_sep, plot.by = "Mac", split.by = "Sample")

dat <- as.data.frame(mac_sep@MacPolarizeR$SingleCell)
dat[,"Sample"] <- mac_sep@MetaData[["Sample"]]
dat[,"Cluster"] <- mac_sep@MetaData[["Cluster"]]

ggplot(dat, aes(x=dat[,1], y=dat[,2])) +
  labs(x="M1", y="M2", color="", fill="") +
  geom_vline(xintercept = 1, colour = 'grey30') +
  geom_hline(yintercept = 1, colour = 'grey30') +
  theme(text = element_text(size = 17),
        panel.background = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.line = element_line(),
        strip.background = element_rect(color = "black")) +
  stat_density_2d(aes(color = dat[,3], fill = dat[,3]), geom = "polygon", alpha = 0.2, position = "identity") +
  #geom_point(aes(color = Sample, fill = Sample), size=3) +
  scale_x_continuous(limits = c(-0.1, 2.5)) +
  scale_y_continuous(limits = c(0, 2.5))

