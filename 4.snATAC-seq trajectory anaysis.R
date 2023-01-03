##snATAC-seq trajectory anaysis.
library(Seurat)
library(Signac) 
library(monocle3)
library(cicero)
library(BuenColors)
library(EnsDb.Mmusculus.v79)
library(here)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
set.seed(1234)

atac<- readRDS("myonuclei_atac.drs")
peaks <- data.frame(atac@assays$peaks@counts@Dimnames[1])
colnames(peaks) <- "coord"
peaks.bed <- data.frame(str_split(peaks$coord, pattern=":|-", simplify=TRUE))
colnames(peaks.bed) <- c("chr","start","end")
write.table(peaks.bed, file = "peaks.bed", row.names = FALSE, quote = FALSE, col.names = FALSE)

# create cell data set object with cicero constructor
input_cds <- as.cell_data_set(atac)
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, num_dim =50)
input_cds <- cluster_cells(input_cds,resolution = 0.0001)
plot_cells(input_cds,color_cells_by = "celltype")
input_cds <- learn_graph(input_cds)
input_cds <- order_cells(input_cds) 
cds_subset <- input_cds[,colData(input_cds)$celltype8 %in% c("TypeIIslow","TypeIIfast")] 
cds_subset <- choose_cells(cds_subset) 
cds_subset <- cluster_cells(cds_subset,resolution = 0.00012)
cds_subset <- learn_graph(cds_subset)
cds_subset <- order_cells(cds_subset) 
plot_cells(cds_subset, color_cells_by = "pseudotime") +  ggtitle("Pseudotime orderings of TypeII snATAC dataset")

cds_subset_lin <- cds_subset[,is.finite(pseudotime(cds_subset))]
pData(cds_subset_lin)$Pseudotime <- pseudotime(cds_subset_lin)
pData(cds_subset_lin)$cell_subtype <- cut(pseudotime(cds_subset_lin), 20)
binned_cds_subset_lin <- aggregate_by_cell_bin(cds_subset_lin, "cell_subtype")

# fit a model to find regions that differ relative to pseudotime
set.seed(1000)
acc_fits <- fit_models(binned_cds_subset_lin[sample(1:nrow(fData(binned_cds_subset_lin)), 206818),], model_formula_str = "~Pseudotime + num_genes_expressed" )
fit_coefs <- coefficient_table(acc_fits)

# Subset out the differentially accessible sites with respect to Pseudotime
pseudotime_terms <- subset(fit_coefs, term == "Pseudotime" & p_value < 0.001)
saveRDS(pseudotime_terms,file = "pseudotime_terms_p0.001.rds")
pseudotime_terms <-pseudotime_terms readRDS("pseudotime_terms.rds")
head(pseudotime_terms)
df <- as.data.frame(pseudotime_terms)
sapply(df, class)
pseudo.peak <- df$site_name
write.csv(pseudo.peak ,file="pseudo.peak.csv")

# Generate a heatmap
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)
library(Matrix)
library(nabor)
library(dplyr)
library(viridis)
library(data.table)
library(scales)
library(mgcv)
library(gplots)
library(RColorBrewer)
library(parallel)
library(speedglm)
library(gtools)
library(uwot)
library(splines)
library(RANN)

pseudo.peak <- read.csv("pseudo.peak.csv")
genes<-pseudo.peak
pt.matrix <-  as.matrix(atac@assays$peaks[match(genes,rownames(atac)),order(atac$pseudo)])
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes;

cols <- viridis(100)
#cols <- colorRampPalette(c("dodgerblue4", "deepskyblue","grey85", "darkorange", "firebrick3"))(100)
#cols <- colorRampPalette(c("grey80", "grey75",brewer.pal(7, "YlGnBu")[2:7]))(100)
cols <- colorRampPalette(c("paleturquoise4", "white","palevioletred3"))(100)
row.o <- apply( pt.matrix, 1, which.max)
fit <- pt.matrix[order(row.o, decreasing=F),]
heatmap.2(fit, trace="none", col=cols, Colv=NA, Rowv=NA, dendrogram="none",
          scale="none", labRow = NA, labCol=NA, useRaster=T, 
          ylab=paste("ACRs", paste0("(n=",nrow(fit),")"), sep=" "))
write.table(fit, file="6_TypeII_pseudotime_ACRs_pt.txt", quote=F, row.names=T, col.names=T, sep="\t")
saveRDS(fit,file="6_TypeII_pseudotime_ACRs_pt.rds")
