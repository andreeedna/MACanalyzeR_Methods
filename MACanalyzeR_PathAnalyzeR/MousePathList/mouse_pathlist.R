setwd("MACanalyzeR_PathAnalyzeR/MousePathList/")

gmtPathways <- function(gmt.file) {
  pathwayLines <- strsplit(readLines(gmt.file), "\t")
  pathways <- lapply(pathwayLines, tail, -2)
  names(pathways) <- sapply(pathwayLines, head, 1)
  pathways
}

## go
mouse_go <- gmtPathways("GO.v5.2.symbols_mouse.gmt")
save(mouse_go, file = "mouse_go.rda")

## kegg
mouse_kegg <- gmtPathways("kegg.v5.2.symbols_mouse.gmt")
save(mouse_kegg, file = "mouse_kegg.rda")

## reactome
mouse_react <- gmtPathways("reactome.v5.2.symbols_mouse.gmt")
save(mouse_react, file = "mouse_reactome.rda")
