library(Seurat)
library(tidyr)

setwd("~/GitHub/MACanalyzeR_MaterialsMethods/FoamSpotteR/FoamSpotteR_scAnalysis/")

#######################################
##  load data
#######################################
## GSE97310 - HEALTY HFD
for(x in list.dirs(path='GSE97310', recursive = F, full.names = F)){
  name <- gsub('','',x)
  cts <- ReadMtx(mtx = paste0("GSE97310/", x, '/matrix.mtx.gz'),
                 features = paste0("GSE97310/", x,'/genes.tsv.gz'),
                 cells = paste0("GSE97310/", x, '/barcodes.tsv.gz'))
  assign(x, CreateSeuratObject(counts = cts, project = "GSE97310"))
  rm(cts, name, x)
}

## already filtered
HEALTY[['percent.mt']] <- PercentageFeatureSet(HEALTY, pattern = '^mt-')
VlnPlot(HEALTY, c("nCount_RNA", "nFeature_RNA", "percent.mt"))
HFD_11W[['percent.mt']] <- PercentageFeatureSet(HFD_11W, pattern = '^mt-')
VlnPlot(HFD_11W, c("nCount_RNA", "nFeature_RNA", "percent.mt"))
HFD_20W[['percent.mt']] <- PercentageFeatureSet(HFD_20W, pattern = '^mt-')
VlnPlot(HFD_20W, c("nCount_RNA", "nFeature_RNA", "percent.mt"))


## GSE116240 - FOAM LDLR
for(x in list.dirs(path='GSE116240', recursive = F, full.names = F)){
  name <- gsub('','',x)
  cts <- ReadMtx(mtx = paste0("GSE116240/", x, '/matrix.mtx.gz'),
                 features = paste0("GSE116240/", x,'/genes.tsv.gz'),
                 cells = paste0("GSE116240/", x, '/barcodes.tsv.gz'))
  assign(x, CreateSeuratObject(counts = cts, project = "GSE116240"))
  rm(cts, name, x)
}

##  already filtered
FOAM[['percent.mt']] <- PercentageFeatureSet(FOAM, pattern = '^mt-')
VlnPlot(FOAM, c("nCount_RNA", "nFeature_RNA", "percent.mt"))
LDLR[['percent.mt']] <- PercentageFeatureSet(LDLR, pattern = '^mt-')
VlnPlot(LDLR, c("nCount_RNA", "nFeature_RNA", "percent.mt"))

####################################################
##  merge data
####################################################
athero <- merge(x= HEALTY , y= c(HFD_11W, HFD_20W, FOAM, LDLR), 
                add.cell.id = c('HEALTY', 'HFD-11W', 'HFD-20W', 'FOAM', 'LDLR'))

athero$sample <- rownames(athero@meta.data)
athero@meta.data <- separate(athero@meta.data, col='sample', into=c('Sample','Barcode'), sep='_')

rm(HEALTY, HFD_11W, HFD_20W, FOAM, LDLR)

####################################################
##  integrate data
####################################################
athero <- NormalizeData(athero)
athero <- FindVariableFeatures(athero,selection.method = "vst", nfeatures = 2000)
athero <- ScaleData(athero, features = VariableFeatures(athero))
athero <- RunPCA(athero, features = VariableFeatures(athero), verbose=F)
athero <- IntegrateLayers(
  object = athero, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = T)

ElbowPlot(athero, reduction = "pca", 50)

athero <- FindNeighbors(athero, reduction = "integrated.cca", dims = 1:20)
athero <- FindClusters(athero, resolution = 0.2)
athero <- RunUMAP(athero, reduction = "integrated.cca", dims = 1:20)

athero <- JoinLayers(athero)

athero$Sample <- factor(athero$Sample, levels = c('HEALTY', 'HFD-11W', 'HFD-20W', 'LDLR', 'FOAM'))
DimPlot(athero, split.by = "Sample")

DimPlot(athero, label = T)

###############
##  marker annotation
###############
##  macrophages
FeaturePlot(athero, c("Mrc1", "Adgre1", "Cd163", "Lyve1"), order = T)
##  monocytes
FeaturePlot(athero, c("Ly6c2", "Ccr2", "Cx3cr1"), order = T)
##  dendridic 
FeaturePlot(athero, c("Ccr9", "Siglech", "Ifi205", "Naaa", "Napsa", "Cd209a"), order = T)
##  B cell
FeaturePlot(athero, c("Ms4a1", "Cd79a"), order = T)
##  T cell
FeaturePlot(athero, c("Cd4", "Cd8a", "Nkg7"), order = T)
##  Neutrophil
FeaturePlot(athero, c("S100a8", "S100a9"), order = T)
##  Endo
FeaturePlot(athero, c("Cav1", "Kdr"), order = T)
##  Endo
FeaturePlot(athero, c("Acta2", "Col1a1"), order = T)

new.id <- c("Macrophages", "Macrophages", "Macrophages", "NK Cells", "Macrophages", "cDC2", "Fibro",
            "T Cells", "B Cells", "cDC1", "T Cells", "Neutrophils", "Endo")
athero$Cells <- factor(athero$seurat_clusters, levels = 0:12, labels = new.id)
Idents(athero) <- "Cells"
rm(new.id)

DimPlot(athero, cols = RColorBrewer::brewer.pal("Set1", n=9))
DimPlot(athero, split.by = "orig.ident", cols = RColorBrewer::brewer.pal("Set1", n=9))
DimPlot(athero, split.by = "Sample", cols = RColorBrewer::brewer.pal("Set1", n=9))
# saveRDS(athero, "totale_annotato.rds")

DotPlot(athero, rev(c("Mrc1", "Adgre1", "Nkg7", "Klrd1", "Napsa", "Cd209a", "Acta2", "Col1a1", "Cxcr6", "Cd3g", "Ms4a1", 
                  "Cd79a", "Ifi205", "Naaa", "S100a8", "S100a9", "Cav1", "Kdr")), cols = "RdYlBu") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

################################
##  macrophages
################################
mac <- subset(athero, idents = "Macrophages")
mac <- NormalizeData(mac)
mac <- FindVariableFeatures(mac, selection.method = "vst", nfeatures = 2000)
mac <- ScaleData(mac)
mac <- RunPCA(mac, features = VariableFeatures(mac), verbose=F)

mac <- IntegrateLayers(
  object = mac, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = T)

ElbowPlot(mac, reduction = "pca", 50)

mac <- FindNeighbors(mac, reduction = "integrated.cca", dims = 1:20)
mac <- FindClusters(mac, resolution = 0.2)
mac <- RunUMAP(mac, reduction = "integrated.cca", dims = 1:20)

mac <- JoinLayers(mac)

DimPlot(mac, group.by = "seurat_clusters", label = T)
DimPlot(mac, split.by = "Sample")

## markers
FeaturePlot(mac, c("Lpl", "Ctsd", "Spp1", "Cd9", "Trem2", "Fabp5"),
            order = T, combine = T, ncol = 3, max.cutoff = "q95", min.cutoff = "q5") &
  scale_color_viridis_c()

new.id <- c("fMAC-", "fMAC-", "fMAC+", "fMAC-", "fMAC+", "fMAC-", "fMAC-",
            "fMAC-", "fMAC-", "fMAC-", "fMAC-")
mac$Cells <- factor(mac$seurat_clusters, levels = 0:10, labels = new.id)
Idents(athero) <- "Cells"
rm(new.id)

DimPlot(mac)

DimPlot(mac, cols = foam_color, pt.size = 1.5)
DimPlot(mac, split.by = "orig.ident", cols = c("fMAC+"="#FF7F00", "fMAC-"="#377EB8"), pt.size = 1.5)
DimPlot(mac, split.by = "Sample", cols = c("fMAC+"="#FF7F00", "fMAC-"="#377EB8"), pt.size = 1.5)

