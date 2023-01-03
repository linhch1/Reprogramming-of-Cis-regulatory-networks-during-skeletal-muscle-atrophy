## Process snATAC datasets 
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(cowplot)
set.seed(1234)
counts <- Read10X_h5("filtered_tf_bc_matrix.h5")
metadata <- read.csv(
  file = "singlecell.csv",
  header = TRUE,
  row.names = 1)

atac<- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = 'fragments.tsv.gz',
  min.cells = 1)

atac<- CreateSeuratObject(
  counts = atac,
  assay = 'peaks',
  project = 'atac',
  meta.data = metadata)

# Extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# Change to UCSC style since the data was mapped to mm10
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# Add the gene information to the object
Annotation(atac) <- annotations
atac<- NucleosomeSignal(object = atac)
atac<- TSSEnrichment(atac, fast = FALSE)
atac$high.tss <- ifelse(atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(atac, group.by = 'high.tss') + NoLegend()
atac$pct_reads_in_peaks <- atac$peak_region_fragments / atac$passed_filters * 100
atac$blacklist_ratio <- atac$blacklist_region_fragments / atac$peak_region_fragments

VlnPlot(
  object = atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5)

atac<- subset(
  x = atac,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 30000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.025 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

atac<- RunTFIDF(atac)
atac<- FindTopFeatures(atac, min.cutoff = 'q0')
atac<- RunSVD(object = atac)

library(harmony)
atac <- RunHarmony(
  object = atac,
  group.by.vars = 'sti',
  reduction = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE
)
atac <- RunUMAP(atac, dims = 2:40, reduction = 'harmony')

atac<- FindNeighbors(
  object = atac,
  reduction = 'lsi',
  dims = 2:40
)
atac<- FindClusters(
  object = atac,
  algorithm = 3,
  resolution = 2,
  verbose = FALSE
)

DimPlot(object = atac, label = TRUE) + NoLegend()

# Compute gene activities
gene.activities <- GeneActivity(atac)

# Add the gene activity matrix to the Seurat object as a new assay
atac[['RNA']] <- CreateAssayObject(counts = gene.activities)
atac<- NormalizeData(
  object = atac,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(atac$nCount_RNA)
)

# Nuclei type annotation
DefaultAssay(atac) <- 'RNA'

# Load the pre-processed scRNA-seq data
pbmc <- readRDS("skeletal_snRNA_pbmc.rds")
pbmc <- FindVariableFeatures(
  object = pbmc,
  nfeatures = 5000
)

transfer.anchors <- FindTransferAnchors(
  reference = pbmc ,
  query = atac,
  reduction = 'cca',
  dims = 2:50
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc$celltype,
  weight.reduction = atac[['lsi']],
  dims = 2:50
)

atac<- AddMetaData(object = atac, metadata = predicted.labels)
FeaturePlot(
  object = atac,
  features = c("Myh7","Myh2","Myh1","Myh4"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 4
)
# Annotate the nuclei based on the prediction score and lineage specific gene activity. 
Celltype<-read.csv(file="celltype.csv")
atac$celltype<-Celltype$celltype
saveRDS(atac, file="snATAC_all_celltype.rds")
Fig1B<-DimPlot(pbmc, group.by = c("celltype"))
Fig1C<-DimPlot(atac, group.by = c("celltype"))
# Select nuclei from control samples 
atac<-subset(atac, sti=="C")
pbmc<-subset(pbmc,sti=="CTL")
Fig1D<-CoveragePlot(
  object = atac,
  region = c("Myh7","Myh2","Myh1","Myh4","Ache","Tnmd","Pax7","Pdgfra","Mrc1","Flt1","Pdgfrb","Lipe"),
  extend.upstream =2000,
  annotation = EnsDb.Mmusculus.v79
)

# Integrating with scRNA-seq data
genes.use <- VariableFeatures(pbmc)
refdata <- GetAssayData(pbmc, assay = "RNA", slot = "data")[genes.use, ]
coembed <- merge(x = pbmc, y = pbmc.atac)

# Visualize the co-embedding of both datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
Fig1E<-DimPlot(coembed, group.by = c("orig.ident"))

## Identify differentially accessible chromatin regions between celltypes
DefaultAssay(atac) <- 'peaks'
GetMarkers <- function(cluster, seurat_aggregate) {
  print(paste0("Finding DAR for: ",cluster))
  atac <- seurat_aggregate
  dar <- FindMarkers(atac, 
                     ident.1 = cluster,    
                     test.use = 'LR', 
                     latent.vars = "peak_region_fragments",
                     min.pct = 0.1,only.pos = T) 
  cf <- ClosestFeature(atac,rownames(dar))
  return(cbind(dar, gene=cf$gene_name, distance=cf$distance,gene_biotype=cf$gene_biotype))
}

# FindMarkers and write to an xlsx file with default parameters
idents <- levels(atac)
list.cluster.dar <- lapply(idents, function(x) {GetMarkers(x, seurat_aggregate = atac)})
write.xlsx(list.cluster.dar, file = "20220831_CTL_atac_cluster_nuclei_dar.xlsx", sheetName = idents, rowNames = T)

darfile <- "20220831_CTL_atac_cluster_nuclei_dar.xlsx"
idents <- getSheetNames(darfile)
list.dar <- lapply(idents, function(x) {
  df <- read.xlsx(darfile, sheet = x, rowNames = T) %>%
    rownames_to_column(var = "coord") %>%
    dplyr::mutate(celltype = x) # annotate each region with its corresponding celltype
})

# Identify all unique cell-type-specific peaks and filter for logfc > 0.25
all_dar <- bind_rows(list.dar) %>%
  dplyr::filter(avg_log2FC > 0.25) %>%dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::select("coord") %>%
  dplyr::distinct()

dar_aver <- AverageExpression(atac, features = all_dar$coord, assays = "peaks")
Fig2A <- pheatmap::pheatmap(dar_aver[["peaks"]],scale = "row" , cluster_rows = F, cluster_cols = F,
                            show_rownames = FALSE) 

Fig2B <- CoveragePlot(
  object = atac,
  region = c("Trim63","Prkag3","Chodl","Scara5"),
  extend.upstream =2000,
  annotation = EnsDb.Mmusculus.v79
)

Fig2C<-FeaturePlot(pbmc, features = c("Trim63","Prkag3","Chodl","Scara5"))

## DARs location presentation
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)
library(org.Mm.eg.db)

# convert the DAR to GRanges objects to annotate
all_dar.gr <- StringToGRanges(all_dar$coord, sep = c(":","-"))
list.dar.gr <- lapply(seq(list.dar), function(x) {
  df <- list.dar[[x]]
  gr <- StringToGRanges(df$coord, sep = c(":","-"))
  return(gr)
})
names(list.dar.gr) <- idents

# annotate the list of GRanges DAR for each cell type
list.peakAnno <- lapply(list.dar.gr, annotatePeak, TxDb = txdb,
                        tssRegion = c(-3000, 3000), verbose = FALSE)
Fig S2A <- plotDistToTSS(list.peakAnno)
Fig S2B <- plotAnnoBar(list.peakAnno) 

# GO enrichment analysis
genes = lapply(list.peakAnno, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compGO<- compareCluster(geneCluster   = genes,
                        fun           = "enrichGO",
                        pvalueCutoff  = 0.05,
                        pAdjustMethod = "BH",OrgDb =org.Mm.eg.db)

FigS2E<-dotplot(compGO, showCategory =5, title = "GO Enrichment Analysis")+RotatedAxis()+ scale_color_viridis(direction = -1)