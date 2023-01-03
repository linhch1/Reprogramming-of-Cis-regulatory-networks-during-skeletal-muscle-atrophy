## Cis-coaccessibility networks generation with cicero
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(cicero)
# Convert to CellDataSet format and make the cicero object
bone.cds <- as.cell_data_set(x = atac)
bone.cicero <- make_cicero_cds(bone.cds, reduced_coordinates = reducedDims(bone.cds)$UMAP)
# Get the chromosome sizes from the Seurat object
genome <- seqlengths(atac)

# Convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# Run cicero
conns <- run_cicero(bone.cicero, genomic_coords = genome.df, sample_num = 100)
ccans <- generate_ccans(conns)
head(ccans)

links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(atac) <- links
saveRDS(atac,file = "atac_with_links.rds")

# Create a column that identifies which connections belong to a CCAN
ccan1 <- left_join(conns, ccans, by=c("Peak1" = "Peak"), all.x=TRUE)
colnames(ccan1)[4] <- "CCAN1"
ccan2 <- left_join(conns, ccans, by=c("Peak2" = "Peak"), all.x=TRUE)
colnames(ccan2)[4] <- "CCAN2"
df <- cbind(ccan1, CCAN2=ccan2$CCAN2) %>%
  dplyr::mutate(CCAN = ifelse(CCAN1 == CCAN2, CCAN1, 0)) %>%
  dplyr::select(-CCAN1, -CCAN2)
fwrite(df, file = "ciceroConns.allcells.csv", row.names = TRUE)

Fig2H(Myh7)<-CoveragePlot(atac,region = c("chr14-54992000-55008000"))#For example

# DARs are linked by the connections and create circos plots
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(openxlsx)
library(circlize)
library(dplyr)
library(tibble)
library(data.table)
library(Signac)
library(stringr)
library(RColorBrewer)

# read in cell-type-specific CCAN and filter by coaccess threshold > 0.2
dar_files <- list.files("ccans", recursive = TRUE, pattern = "ciceroConns", full.names = TRUE)
list.dar <- lapply(dar_files, function(file_path) {
  fread(file_path) %>%
    dplyr::filter(coaccess > 0.2) %>%
    dplyr::select(Peak1, Peak2)
})


idents<-c("Adipocyte","Endothelial","FAPs","Macrophage","MTJ","MuSCs","Pericyte","TypeIIfast","TypeIIslow")
# convert the DAR to GRanges objects to annotate
list.dar.peak1.gr <- lapply(seq(list.dar), function(x) {
  df <- list.dar[[x]]
  gr <- StringToGRanges(df$Peak1, sep = c("-","-"))
  return(gr)
})
names(list.dar.peak1.gr) <- idents

list.dar.peak2.gr <- lapply(seq(list.dar), function(x) {
  df <- list.dar[[x]]
  gr <- StringToGRanges(df$Peak2, sep = c("-","-"))
  return(gr)
})
names(list.dar.peak2.gr) <- idents

# annotate the list of GRanges DAR for each cell type with the peak location
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
list.peak1.Anno <- lapply(list.dar.peak1.gr, annotatePeak, TxDb = txdb,
                          tssRegion = c(-3000, 3000), verbose = FALSE)
list.peak1.loc <- lapply(seq(list.peak1.Anno), function(x) {
  loc <- list.peak1.Anno[[x]]@anno$annotation 
  loc <- str_split(loc, pattern = " ", simplify = TRUE)[,1]
  loc <- str_replace(loc, pattern="3'", replacement="3p_UTR")
  loc <- str_replace(loc, pattern="5'", replacement="5p_UTR")
  loc <- str_replace(loc, pattern="Distal", replacement="Intergenic")
  return(loc)
})

list.peak2.Anno <- lapply(list.dar.peak2.gr, annotatePeak, TxDb = txdb,
                          tssRegion = c(-3000, 3000), verbose = FALSE)
list.peak2.loc <- lapply(seq(list.peak2.Anno), function(x) {
  loc <- list.peak2.Anno[[x]]@anno$annotation 
  loc <- str_split(loc, pattern = " ", simplify = TRUE)[,1]
  loc <- str_replace(loc, pattern="3'", replacement="3p_UTR")
  loc <- str_replace(loc, pattern="5'", replacement="5p_UTR")
  loc <- str_replace(loc, pattern="Distal", replacement="Intergenic")
  return(loc)
})

# create df with peak1 and peak2 locations
list.peaks.loc.df <- lapply(seq(idents), function(x) {
  df <- cbind(list.peak1.loc[[x]], list.peak2.loc[[x]]) %>%
    as.data.frame()
  colnames(df) <- c("Peak1","Peak2")
  counts <- dplyr::count(df, Peak1, Peak2) %>% # count the number of times each pair occurs
    as.data.frame()
})
names(list.peaks.loc.df) <- idents

idents
# generate circos plots of predicted chromatin-chromatin interactions and save to plots directory
Plot_Cicero_Anno <- function(ident, list.peaks.loc.df) {
  clusterID <- ident
  toplot <- as.data.frame(list.peaks.loc.df[[ident]])
  colnames(toplot) <- c("Peak1","Peak2","n")
  print(clusterID)
  # convert to an adjacency list with a value indicating the number of connections between
  # each of the unique genomic location pairs
  unique_combos <- !duplicated(t(apply(toplot, 1, sort)))
  toplot <- toplot[unique_combos, ]
  toplot <- toplot[order(toplot$n), ]
  # toplot$n <- log2(toplot$n)
  
  # grid.col = c('3p_UTR'="red", '5p_UTR'="black", Distal_Intergenic="blue",
  # Downstream="grey", Exon="purple", Intron="orange", Promoter="green")
  grid.col = brewer.pal(7, "Paired")
  
  circos.clear()
  pdf(paste0("plots/ccans.annotation.",clusterID,".pdf"), width=10, height=5)
  par(mar=c(0.5,0.5,0.5,0.5))
  circos.par(gap.after = 5)
  chordDiagram(toplot,
               grid.col = grid.col,
               annotationTrack = "grid",
               preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(toplot))))))
  
  # we go back to the first track and customize sector labels
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter,
                CELL_META$ylim[1],
                CELL_META$sector.index,
                facing = "clockwise",
                niceFacing = TRUE,
                adj = c(0, 0.5))
  }, bg.border = NA) # here set bg.border to NA is important
  dev.off()
}

# save circlize plots to file
lapply(idents, function(ident) {
  Plot_Cicero_Anno(ident, list.peaks.loc.df=list.peaks.loc.df)
})