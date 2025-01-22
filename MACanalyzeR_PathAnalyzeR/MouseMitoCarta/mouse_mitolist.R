## creation of PathAnalyzeR compatible list starting from MitoCarta xlsx
setwd("MACanalyzeR_PathAnalyzeR/MouseMitoCarta/")

mouse_MitoCarta <- list()

for (i in 1:8) {
  mitopath <- openxlsx::read.xlsx("MouseMitoCartaR.xlsx", sheet = i, colNames = F, rowNames = F)
  mitolist <- list()
  
  for (m in 1:nrow(mitopath)) {
    mitolist[mitopath[m,1]] <- as.vector(strsplit(mitopath[m,2], split = ", "))
  }

  if (i==1) {
    x <- names(mitolist)
    mouse_MitoCarta[["Total"]] <- mitolist
  } else {
    mouse_MitoCarta[[x[i-1]]] <- mitolist 
  }
}

names(mouse_MitoCarta)

save(mouse_MitoCarta, file = "mouse_MitoCarta.rda")
