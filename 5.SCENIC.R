## SCENIC
library(foreach)
library(SCENIC)
library(Seurat)
library(GENIE3)
library(AUCell)
library(RcisTarget)
library(loomR)

dbFiles <- c(''https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather'', ''https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather'')

download.file(dbFiles[1],'SCENIC/cisTarget_databases/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather')
download.file(dbFiles[2],'SCENIC/cisTarget_databases/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather')
exprMat <- pbmc@assays$RNA@counts
exprMat <- as.matrix(exprMat)
cellInfo <- data.frame(pt@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=='orig.ident')] <- 'sample'
colnames(cellInfo)[which(colnames(cellInfo)=='seurat_clusters')] <- 'cluster'
colnames(cellInfo)[which(colnames(cellInfo)=='sub_type')] <- 'celltype'
cellInfo <- cellInfo[,c('sample','cluster','celltype','Group')]
saveRDS(cellInfo, file='int/cellInfo.Rds')

org="mgi" # or hgnc, or dmel
dbDir="cisTarget_databases" # RcisTarget databases location
myDatasetTitle="SCENIC example" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)

# Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

# Gene regulatory network 
#  (Adjust minimum values according to your dataset)
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)

runGenie3(exprMat_filtered_log, scenicOptions)

# Build and score the GRN 
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] 
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions) 
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

# Binarize network activity 
runSCENIC_4_aucell_binarize(scenicOptions)

# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/):
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))
par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)

# Export:
export2scope(scenicOptions, exprMat)

# Exploring output
# Check files in folder 'output'
# .loom file @ http://scope.aertslab.org
# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Sox8"]
viewMotifs(tableSubset)
# output/Step2_regulonTargetsInfo.tsv in detail:
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Stat6" & highConfAnnot==TRUE]
viewMotifs(tableSubset)