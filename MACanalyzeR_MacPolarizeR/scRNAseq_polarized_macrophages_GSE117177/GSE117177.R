## scRNAseq analysis of Bone Marrow Macrophages after M1 and M2 stimuli
## analysis performed with Seurat 4
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)

setwd("MACanalyzeR_MacPolarizeR/scRNAseq_polarized_macrophages_GSE117177/")

######################
# LOAD DATA
######################
for(x in list.dirs(path='.', recursive = F, full.names = F)){
  name <- gsub('','',x)
  cts <- ReadMtx(mtx = paste0(x, '/matrix.mtx.gz'),
                 features = paste0(x, '/genes_adj.tsv.gz'),
                 cells = paste0(x, '/barcodes.tsv.gz'))
  assign(x, CreateSeuratObject(counts = cts, project = x))
  rm(cts, name, x)
}

###############################
# QC and FILTERING
###############################
#calculate % MT genes
M0[['percent.mt']] <- PercentageFeatureSet(M0, pattern = '^[Mm][Tt]-')
M1[['percent.mt']] <- PercentageFeatureSet(M1, pattern = '^[Mm][Tt]-')
M2[['percent.mt']] <- PercentageFeatureSet(M2, pattern = '^[Mm][Tt]-')

# M0
M0_filt <- subset(M0, subset= percent.mt < 5)
VlnPlot(M0_filt, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
M0_filt

# M1
M1_filt <- subset(M1, subset= percent.mt < 5)
VlnPlot(M1_filt, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
M1_filt

# M2
M2_filt <- subset(M2, subset= percent.mt < 5)
VlnPlot(M2_filt, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
M2_filt


#sample the same number of cells - campione minore WT16w
set.seed(25)
n.cell.sample <- 4000

m0.sampled.cells <- sample(rownames(M0_filt@meta.data), size = n.cell.sample, replace = F)
M0_filt <- subset(M0_filt, cells = m0.sampled.cells)
M0_filt
rm(m0.sampled.cells)

m1.sampled.cells <- sample(rownames(M1_filt@meta.data), size = n.cell.sample, replace = F)
M1_filt <- subset(M1_filt, cells = m1.sampled.cells)
M1_filt
rm(m1.sampled.cells)

m2.sampled.cells <- sample(rownames(M2_filt@meta.data), size = n.cell.sample, replace = F)
M2_filt <- subset(M2_filt, cells = m2.sampled.cells)
M2_filt
rm(m2.sampled.cells)
rm(n.cell.sample)

rm(M0, M1, M2)


#############################
# MERGE SEURAT OBJECTS
#############################
bey <- merge(x= M0_filt , y=c(M1_filt, M2_filt),
                       add.cell.id = c('M0', 'M1', 'M2'))

bey$sample <- rownames(bey@meta.data)
bey@meta.data <- separate(bey@meta.data, col='sample',
                                    into=c('Sample','Barcode'),
                                    sep='_')

rm(M0_filt, M1_filt, M2_filt)

############################
# NON INTEGRATION STEP
############################
bey <- NormalizeData(bey)
bey <- FindVariableFeatures(bey, selection.method = "vst", nfeatures = 2000)
bey <- ScaleData(bey, features = rownames(bey))
bey <- RunPCA(bey, features = VariableFeatures(object = bey))

ElbowPlot(bey, 50)

bey <- RunTSNE(object = bey, reduction = 'pca', dims= 1:10)
bey <- FindNeighbors(bey, dims = 1:10)
bey <- FindClusters(bey, resolution = 0.2)

#visualize
DimPlot(bey, reduction = 'tsne', group.by = 'Sample')
DimPlot(bey, reduction = 'tsne', label = TRUE)


#############################
##  MACanalyzeR
#############################
library(MACanalyzeR)
mac_bey <- CreateMacObj(bey, "Sample", "Sample", "tsne")
mac_bey <- MacPolarizeR(mac_bey)
MacPlot(mac_bey, plot.by = "Mac") + MacPlot(mac_bey)

dat <- as.data.frame(mac_bey@MacPolarizeR$SingleCell)
dat[,"Sample"] <- mac_bey@MetaData[["Sample"]]
dat[,"Cluster"] <- mac_bey@MetaData[["Cluster"]]


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

