library(Seurat)
library(swne)
library(cellMapper)

## Set working directory to file directory
setwd("teratoma-analysis-code/Figure1")

## Input/output files
input.file <- "dm-ter-mouse-merged_clustering_v3.seurat.Robj"
output.file <- gsub("clustering_v3.seurat.Robj", "mapping_v3.RData", input.file)

## These reference files were generated from the Mouse Cell Atlas data
## http://bis.zju.edu.cn/MCA/
## These paths may change depending on where you download your data
ref.seurat.file <- "../Reference_Data/mca_seurat_mouse.RData"
ref.markers.file <- "../Reference_Data/mca_markers_mouse.txt"


#### Input teratoma clustering Seurat ####

## Load reference data
load(ref.seurat.file)

## Load teratoma data
ter <- readRDS(input.file)

## Select genes to train on
ref.markers.df <- read.table(ref.markers.file, header = T, sep = " ")
ter.markers.df <- read.table(gsub(".seurat.Robj", ".markers.tsv", input.file), header = T, sep = "\t")

ref.norm.counts <- ScaleCounts(mca@raw.data[,mca@cell.names])
ter.norm.counts <- GetAssayData(ter, assay = "integrated", slot = "scale.data")

genes.use <- unique(c(as.character(ref.markers.df$gene), VariableFeatures(ter)))
genes.use <- Reduce(intersect, list(genes.use, rownames(ref.norm.counts), rownames(ter.norm.counts)))
length(genes.use)

## Run PCA on reference data
pcs.use <- 50
ref.ident <- as.character(mca@ident); names(ref.ident) <- mca@cell.names;
ref.ident <- sapply(ref.ident, ExtractField, 1, delim = "_")
ref.ident <- sapply(ref.ident, ExtractField, 1, delim = "\\(")
ref.ident <- factor(ref.ident)
table(ref.ident)

train.knn <- TrainKNN(ref.norm.counts, genes.use, ref.ident, pcs.use = pcs.use)

## Classify cells with kNN classifier
cell.mapping <- MapKNN(ter.norm.counts, train.knn, genes.use = genes.use, use.pca = T, k = 50)

## Summarize according to teratoma cluster
ter.ident <- Idents(ter)
cl.mapping <- SummarizeKNN(cell.mapping, ter.ident, ident.thresh = 0.1, assign.cells = T, min.ident.frac = 0.1)
colnames(cl.mapping) <- paste0("C", colnames(cl.mapping))

pdf(gsub("clustering_v3.seurat.Robj", "tissue_mapping_v3.pdf", input.file), width = 9, height = 6)
ggHeat(cl.mapping, clustering = "none", y.lab.size = 10)
dev.off()

## Look at top markers
top.markers.df <- subset(ter.markers.df, avg_logFC > 0.25)
top.markers.df <- Reduce(rbind, by(top.markers.df, top.markers.df$cluster, head, n = 10))

## Create PDF plots
umap.emb <- Embeddings(ter, reduction = "umap")

cluster.mapping <- c("0"="Stromal", "1"="Stromal", "6"="Stromal", "9"="Stromal",
                     "11"="Stromal", "2"="Macrophage", "3"="Macrophage", "4"="Macrophage",
                     "8"="Macrophage", "12"="Macrophage", "16"="Macrophage", "5"="Endothelial",
                     "7"="Dendritic", "10"="Dendritic", "13"="Smooth Muscle",
                     "14"=NA, "15"=NA, "17"=NA)
named.ident <- plyr::revalue(ter.ident, replace = cluster.mapping)
table(named.ident)

pdf("dm-ter-mouse-merged_umap_named_clusters.pdf", width = 6, height = 6)
PlotDims(umap.emb, sample.groups = named.ident, show.axes = F, show.legend = F,
         pt.size = 0.4, alpha = 0.5, seed = 12535, label.size = 5)
dev.off()

ter.id <- factor(ter$batch)
pdf("dm-ter-mouse-merged_umap_teratoma.pdf", width = 6.5, height = 6)
PlotDims(umap.emb, sample.groups = ter.id, show.axes = F, show.legend = T,
         pt.size = 0.1, alpha = 0.5, seed = 12535, do.label = F, use.brewer.pal = F)
dev.off()

barcoding <- factor(ter$barcoded)
pdf("dm-ter-mouse-merged_umap_barcoding.pdf", width = 6.5, height = 6)
PlotDims(umap.emb, sample.groups = barcoding, show.axes = F, show.legend = T,
         pt.size = 0.1, alpha = 0.5, seed = 12535, do.label = F, use.brewer.pal = F)
dev.off()

gene <- "PAX6"
gene.expr <- ter.norm.counts[gene,]
# pdf(paste0("dm-ter-mouse-merged_umap_", gene, "expr.pdf"), width = 6.5, height = 6)
FeaturePlotDims(umap.emb, gene.expr, feature.name = gene, show.axes = F, pt.size = 0.5,
                alpha.plot = 0.5)
# dev.off()

## Save results
save.image(output.file)
