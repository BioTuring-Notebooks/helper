EnsembleToSymbolMouse <- function(obj) {
    gene.ref <- ensembldb::select(EnsDb.Mmusculus.v79, keys=as.character(rownames(obj)), 
                              keytype = "GENEID", columns = c("SYMBOL","GENEID"))
    genesymbol.ref <- rownames(clp)
    genesymbol.ref[match(gene.ref[[2]], genesymbol.ref)] <- gene.ref[[1]]
    genesymbol.ref <- as.character(genesymbol.ref)
    obj <- RenameGenesSeurat(obj, genesymbol.ref)
    return(obj)
}

RenameGenesSeurat <- function(obj, newnames) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  RNA <- obj@assays$RNA

  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}

PlotCRADGrid <- function(loc) {
    require(ggplot2)
    p <- ggplot(loc, 
       aes(x = x, y = y)) + geom_point(shape=22,color = "#7dc7f5")+
theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    legend.position="bottom",
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5))
    return(p)
}
