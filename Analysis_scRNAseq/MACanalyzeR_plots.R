library(Seurat)
library(MACanalyzeR)
library(dittoSeq)
library(ggplot2)

#######################
##  TOTAL
#######################
## set palette
tot_color_sample <-c("wt"="#E69F00", "hfd"="#56B4E9", "db/db"="#009E73")

adipo_pal <- scales::hue_pal()(5)
names(adipo_pal) <- c("FAP", "Macrophages", "Endothelial Cells", "ASC", "T Cells")

##  dimplot
SCpubr::do_DimPlot(adipo, legend.position = "right")
SCpubr::do_DimPlot(adipo, legend.position = "right", split.by = "Sample", group.by = "Cells", ncol = 4)

##  barplot
dittoBarPlot(adipo, "Sample", "Cells", retain.factor.levels = T)
dittoBarPlot(adipo, "Cells", "Sample", retain.factor.levels = T, color.panel = adipo_pal)

##  heatmap
totmarker <- c("Lum", "Podn", "Col5a3", #FAP,
               "Adgre1", "Itgam", "Lyz2", #mono/MACs
               'Cldn5', 'Cdh5', 'Pecam1', #endothelials
               "Dpp4", "Sema3c", "Scara5", #ASC
               'Trbc2', 'Cd3d', 'Cd3g' #Tcells
)

markexpr <- apply(adipo@assays$RNA$counts[totmarker,], 1, function(x)by(x, adipo$Cells,mean))
pheatmap::pheatmap(t(markexpr), scale = "row", 
                   cluster_rows = F, cluster_cols = F, 
                   color = colorRampPalette(c("blue","white", "red"))(50),
                   border_color = F)


################################
## MACROPHAGES
################################
mac_color <- c("Monocytes"="#3468C0", "PVM"="#86A7FC", "LAM"="#FF9843", "P_LAM"="#FFDD95")

##  dimplot
SCpubr::do_DimPlot(mac, colors.use = mac_color, legend.position = "right")
SCpubr::do_DimPlot(mac, colors.use = mac_color, legend.position = "right", split.by = "Sample", group.by = "Cells", ncol = 4)

##  barplot
dittoSeq::dittoBarPlot(mac, "Cells", "Sample", retain.factor.levels = T, color.panel = mac_color)

##  Marker Heatmap
sub <- c("Mrc1", "Lyve1", "Cd163" #PVM
         ,"Cd36", "Trem2", "Cd9" #LAM 
         ,"Pola1", "Kif11", "Kif15", #P-LAM
         'Ccr2', 'Cd44', 'Cx3cr1' #mono
)

markexpr <- apply(mac@assays$RNA$counts[sub,], 1, function(x)by(x, mac$Cells,mean))
markexpr <- markexpr[c("PVM", "LAM", "P_LAM", "Monocytes"),]
pheatmap::pheatmap(t(markexpr), scale = "row", 
                   cluster_rows = F, cluster_cols = F, 
                   color = colorRampPalette(c("blue","white", "red"))(50),
                   border_color = F)


## Lipid/Lyso/Mito Heatmap
lyso_gen <- matrix(unlist(strsplit(read.table("KEGG_lysosome.txt", header = F, sep = "\t")[,2], ";")), ncol = 2, byrow = T)[,1]
lyso_gen <- intersect(lyso_gen, rownames(mac@assays$RNA$counts))

oxphos_gen <- matrix(unlist(strsplit(read.table("KEGG_OXPHOS.txt", header = F, sep = "\t")[,2], ";")), ncol = 2, byrow = T)[,1]
oxphos_gen <- intersect(oxphos_gen, rownames(mac@assays$RNA$counts))

lipid_gen <-  matrix(unlist(strsplit(read.table("KEGG_LipidMetabolism.txt", header = F, sep = "\t")[,2], ";")), ncol = 2, byrow = T)[,1]
lipid_gen <- intersect(lipid_gen, rownames(mac@assays$RNA$counts))

lyso_gen <- setdiff(setdiff(lyso_gen, lipid_gen), oxphos_gen) 
lipid_gen <- setdiff(setdiff(lipid_gen, lyso_gen), oxphos_gen) 
oxphos_gen <- setdiff(setdiff(oxphos_gen, lipid_gen), lyso_gen) 

gen <- c(lipid_gen, lyso_gen, oxphos_gen)

t <- c(length(intersect(colnames(mexpr), lipid_gen)), 
       length(intersect(colnames(mexpr), lyso_gen)), 
       length(intersect(colnames(mexpr), oxphos_gen)))

# Sample Heatmap
mexpr <- apply(mac@assays$RNA$counts[gen,], 1, function(x)by(x, mac$Sample,mean))
mexpr <- mexpr[,which(colSums(mexpr)!=0)]

ann_row <- data.frame("Sample" = rep(c("wt", "hfd", "db/db"), times =c(1,1,1)))
rownames(ann_row) <- rownames(mexpr)

ann_col <- data.frame("Pathway" = rep(c("Lipid", "Lysosome", "OXPHOS"), times = t))
rownames(ann_col) <- colnames(mexpr)

ann_color <- list("Sample" = c(tot_color_sample[1], tot_color_sample[2], tot_color_sample[3]),
                  "Pathway" = c("Lipid"="gold", "Lysosome"="darkgreen", "OXPHOS"="orange"))

pheatmap::pheatmap(mexpr, scale = "column", 
                   annotation_col = ann_col, annotation_names_col = F,
                   annotation_row = ann_row, annotation_names_row = F,
                   annotation_colors = ann_color, annotation_legend = T,
                   cluster_rows = F, cluster_cols = F, show_colnames = F, 
                   color = viridis::viridis(100), show_rownames = F,
                   border_color = "none", fontsize = 15)

# Cells
mexpr <- apply(mac@assays$RNA$counts[gen,], 1, function(x)by(x, mac$Cells,mean))
mexpr <- mexpr[,which(colSums(mexpr)!=0)]

ann_row <- data.frame("Cells" = rownames(mexpr))
rownames(ann_row) <- rownames(mexpr)

ann_col <- data.frame("Pathway" = rep(c("Lipid", "Lysosome", "OXPHOS"), times = t))
rownames(ann_col) <- colnames(mexpr)

ann_color <- list("Cells" = c("PVM"= "#86A7FC", "LAM"= "#FF9843", "Monocytes"="#3468C0", "P_LAM"= "#FFDD95"),
                  "Pathway" = c("Lipid"="gold", "Lysosome"="darkgreen", "OXPHOS"="orange"))

pheatmap::pheatmap(mexpr, scale = "column", 
                   annotation_col = ann_col, annotation_names_col = F,
                   annotation_row = ann_row, annotation_names_row = F,
                   annotation_colors = ann_color, annotation_legend = T,
                   cluster_rows = F, cluster_cols = F, show_colnames = F, 
                   color = viridis::viridis(100), show_rownames = F,
                   border_color = "none", fontsize = 15)
