## Motif enrichment analysis.
library(Seurat) 
library(Signac) 
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2018)
library(TFBSTools)
library(chromVAR)
library(openxlsx)
library(here)
library(patchwork)
library(motifmatchr)
set.seed(1234)

# Get a list of motif position weight matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2018,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = StringToGRanges(rownames(atac), sep = c(":", "-")),
  pwm = pfm,
  genome = 'BSgenome.Mmusculus.UCSC.mm10',
  sep = c(":", "-")
)

# Create a new Mofif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)

# Add the Motif object to the assay
atac[['peaks']] <- AddMotifObject(
  object = atac[['peaks']],
  motif.object = motif
)

atac <- RegionStats(
  object = atac,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  sep = c(":", "-")
)

# compute motif activities using chromvar
atac <- RunChromVAR(
  object = atac,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  motif.matrix = motif.matrix
)
saveRDS(atac, file = "chromVar.3atac_atac.rds")


GetMotifs <- function(cluster, seurat_aggregate) {
  print(paste0("Finding motifs for: ",cluster))
  atac <- seurat_aggregate
  DefaultAssay(atac) <- 'peaks'
  dac <- FindMarkers(atac,
                     ident.1 = cluster,
                     test.use = 'LR',
                     latent.vars = "nCount_peaks") # find all cluster-specific degs
  enriched.motifs <- FindMotifs(object = atac, features = rownames(dac[dac$p_val < 0.05, ]))
  return(enriched.motifs)
}

GetChromvarActivities <- function(cluster, seurat_aggregate, motif) {
  print(paste0("Finding chromVAR activities for: ",cluster))
  atac <- seurat_aggregate
  DefaultAssay(atac) <- 'chromvar'
  dam <- FindMarkers(atac, 
                     ident.1 = cluster,
                     test.use = 'LR',
                     latent.vars = "nCount_peaks",
                     logfc.threshold = 0,) # find all cluster-specific degs
  motifLookup <- rownames(dam)
  motifNames <- sapply(motifLookup, function(x) motif@motif.names[[x]])
  return(cbind(dam, gene = motifNames))
}

# FindMarkers and write to an xlsx file with default parameters
idents <- levels(atac)
list.cluster.dac <- lapply(idents, function(x) GetMotifs(x, seurat_aggregate = atac))
write.xlsx(list.cluster.dac, file = "20220908_atac_cluster_motifs.xlsx", sheetName = idents, rowNames = T)

list.cluster.dam <- lapply(idents, function(x) GetChromvarActivities(x, seurat_aggregate = atac, motif = motif))
write.xlsx(list.cluster.dam, file = "20220908_atac_cluster_motifs_chromVAR.xlsx", sheetName = idents, rowNames = T)

# Find nuclei type specific TF and make a unique list
tf_celltype.df <- FindAllMarkers(atac, only.pos = T,logfc.threshold = 0.25, assay = "chromvar")
Selected_motifs<-c("MA0838.1","MA0102.4","MA1637.1","MA0154.4","MA1521.1","MA1120.1","MA0687.1","MA0598.3","MA1476.1","MA0841.1","MA1641.1","MA0499.2","MA1525.1","MA0625.1","MA0783.1","MA0796.1")
aver_chromvar <- AverageExpression(atac, assays = "chromvar", features = Selected_motifs)
aver_chromvar <- aver_chromvar[do.call(order, c(aver_chromvar, list(decreasing=TRUE))),]
aver_chromvar$max <- max.col(aver_chromvar)
aver_chromvar <- aver_chromvar[order(aver_chromvar$max), ]
aver_chromvar <- dplyr::select(aver_chromvar, -max)
colnames(aver_chromvar) <- levels(Idents(atac))
library(pheatmap)
paletteLength <- 50
myColor <- viridis::viridis(paletteLength)
Fig2E<-pheatmap(aver_chromvar,scale = "row", cluster_cols=F,cluster_rows = F,show_rownames=F,color = myColor)



Fig2F <- Footprint(
  object = atac,
  motif.name = c("TGIF1", "MYF5", "NFE2","EHF"),
  genome = BSgenome.Mmusculus.UCSC.mm10
)
