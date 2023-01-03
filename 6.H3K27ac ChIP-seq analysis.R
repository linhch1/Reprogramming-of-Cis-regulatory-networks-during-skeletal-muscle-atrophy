## H3K27ac ChIP-seq analysis

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)
library(PeakSegDisk)
peak <- readPeakFile("Sham1_peaks.narrowPeak")
promoter <- getPromoters(TxDb=txdb, upstream=20000, downstream=20000)
newStyle <- mapSeqlevels(seqlevels(peak), "UCSC")
newStyle<-newStyle[!is.na(newStyle)]
peakchr <- renameSeqlevels(peak, newStyle)
tagMatrix <- getTagMatrix(peakchr, windows=promoter)
peak1 <- readPeakFile("deN2_peaks.narrowPeak")
newStyle1 <- mapSeqlevels(seqlevels(peak1), "UCSC")
newStyle1<-newStyle[!is.na(newStyle1)]
peakchr1 <- renameSeqlevels(peak1, newStyle1)
tagMatrix1 <- getTagMatrix(peakchr1, windows=promoter)
files<-list(peak,peak1)
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
tagMatrixList <-list(tagMatrix,tagMatrix1)
peakHeatmap(tagMatrixList,weightCol="V5", upstream=20000, downstream=20000, TxDb=txdb, color="red")
# average profiles of ChIP peaks among different experiments
plotAvgProf(tagMatrixList, xlim=c(-20000,20000))