## creation of PathAnalyzeR compatible list starting from MitoCarta xlsx
setwd("MACanalyzeR_PathAnalyzeR/HumanMitoCarta/")

human_MitoCarta <- list()

for (i in 1:8) {
  mitopath <- openxlsx::read.xlsx("HumanMitoCartaR.xlsx", sheet = i, colNames = F, rowNames = F)
  mitolist <- list()
  
  for (m in 1:nrow(mitopath)) {
    mitolist[mitopath[m,1]] <- as.vector(strsplit(mitopath[m,2], split = ", "))
  }

  if (i==1) {
    x <- names(mitolist)
    human_MitoCarta[["Total"]] <- mitolist
  } else {
    human_MitoCarta[[x[i-1]]] <- mitolist 
  }
}

names(human_MitoCarta[[2]])

save(human_MitoCarta, file = "human_MitoCarta.rda")
