library(tidyverse)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)

setwd("~/GitHub/MACanalyzeR_MaterialsMethods/Analysis_proteomics/")

##################
# LOAD DATA
##################
#LOAD DATA FROM TSV
data <- openxlsx::read.xlsx("BAT_PpargKO_HFD.xlsx")
data <- openxlsx::read.xlsx("BAT_PpargKO_HFD.xlsx", sheet = 2)
head(data)

#CREAdata#CREATE MATRIX 
countdata <- as.matrix(data[,c(9:16,1:8)])
countdata <- as.matrix(data[,c(9:16,5:8)])
head(countdata)
mode(countdata) <- "numeric"
rownames(countdata) <- data$`T:.PG.Genes`
head(countdata)

######################
## HEATMAP PLOTTING
######################
#genes
thermo_genes <- c("Ucp1", "Sdhb", "Sdhc", "Cox7b", "Cox5a", "Cd36", "Idh2")

t <- grep(pattern = "Nduf", rownames(countdata), value = T)
t <- grep(pattern = "Cox", rownames(countdata), value = T)
t <- grep(pattern = "Atp", rownames(countdata), value = T)
t <- grep(pattern = "Slc", rownames(countdata), value = T)

ann_col <- data.frame("Sample" = rep(c("WT", "WT HFD", "PpargKO", "PpargKO HFD"), times =c(4,4,4,4)))
ann_col <- data.frame("Sample" = rep(c("WT", "WT HFD", "PpargKO HFD"), times =c(4,4,4)))
rownames(ann_col) <- colnames(countdata)

ann_color <- list("Sample" = c("PpargKO" = "#E65C19", "PpargKO HFD"="#F8D082", "WT"="#640D6B", "WT HFD"="#B51B75"))

pheatmap(countdata[thermo_genes,], 
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

pheatmap(countdata[t,], 
         cluster_rows = T, 
         cluster_cols = F, 
         scale = "row",
         show_rownames = T, 
         show_colnames = F,
         annotation_col = ann_col,
         annotation_colors = ann_color, 
         annotation_names_row = F,
         annotation_names_col = F, border_color = F,
         color = colorRampPalette(c("blue", "white", "red"))(100))

