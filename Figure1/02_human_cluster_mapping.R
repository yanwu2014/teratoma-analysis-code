library(Seurat)
library(swne)
library(cellMapper)
library(ggplot2)

## Set working directory to file directory
setwd("teratoma-analysis-code/Figure1")

## Input/output files
input.file <- "dm-ter-human-merged_clustering_v3.seurat.Robj"
output.file <- "dm-ter-human-merged_mapping_v3.RData"

## These reference files were generated from the Mouse Cell Atlas data
## http://bis.zju.edu.cn/MCA/
## These paths may change depending on where you download your data to
ref.seurat.file <- "../Reference_Data/mca_seurat_human.RData"
ref.markers.file <- "../Reference_Data/mca_markers_human.txt"


## Load reference data
load(ref.seurat.file)

## Load teratoma data
ter <- readRDS(input.file)

## Select genes to train on
ref.markers.df <- read.table(ref.markers.file, header = T, sep = " ")
ter.markers.df <- read.table(gsub(".seurat.Robj", ".markers.tsv", input.file), header = T, sep = "\t")

ref.norm.counts <- ScaleCounts(mca@raw.data[,mca@cell.names])
ter.norm.counts <- GetAssayData(ter, assay = "integrated", slot = "scale.data")

genes.use <- unique(c(as.character(ref.markers.df$gene), as.character(ter.markers.df$gene), VariableFeatures(ter)))
genes.use <- Reduce(intersect, list(genes.use, rownames(ref.norm.counts), rownames(ter.norm.counts)))
length(genes.use)

## Run PCA on reference data
pcs.use <- 40
# ref.ident <- mca@ident
ref.ident <- as.character(mca@ident); names(ref.ident) <- mca@cell.names;
ref.ident <- sapply(ref.ident, ExtractField, 1, delim = "_")
ref.ident <- sapply(ref.ident, ExtractField, 1, delim = "\\(")
ref.ident <- factor(ref.ident)
table(ref.ident)

train.knn <- TrainKNN(ref.norm.counts, genes.use, ref.ident, pcs.use = pcs.use)

## Classify cells with kNN classifier
cell.mapping <- MapKNN(ter.norm.counts, train.knn, genes.use = genes.use, use.pca = T, k = 40)

## Summarize according to teratoma cluster
ter.ident <- Idents(ter)
cl.mapping <- SummarizeKNN(cell.mapping, ter.ident, ident.thresh = 0.1, assign.cells = T, min.ident.frac = 0.2)
colnames(cl.mapping) <- paste0("C", colnames(cl.mapping))

# pdf(gsub("clustering_v3.seurat.Robj", "tissue_mapping_v3.pdf", input.file), width = 9.5, height = 5)
ggHeat(cl.mapping, clustering = "none", y.lab.size = 10)
# dev.off()

## Clean up some stuff
rm(mca); gc();


## Export top markers for detailed analysis
top.markers.df <- subset(ter.markers.df, avg_logFC > 0.25)
top.markers.df <- Reduce(rbind, by(top.markers.df, top.markers.df$cluster, head, n = 10))
trrust.df <- read.table("trrust_human.tsv", header = F, sep = "\t")
tf.markers.df <- subset(ter.markers.df, gene %in% trrust.df$V1 & avg_logFC > 0.25)
tf.markers.df <- Reduce(rbind, by(tf.markers.df, tf.markers.df$cluster, head, n = 10))

write.table(top.markers.df, file = "dm-ter-human-merged_top_markers.tsv", sep = "\t", row.names = F)
write.table(tf.markers.df, file = "dm-ter-human-merged_TF_markers.tsv", sep = "\t", row.names = F)

## Create PDF plots
umap.emb <- Embeddings(ter, reduction = "umap")

## Seurat cluster to initial cell type mapping. 
## The exact mapping may change due to randomness in the clustering method 
cluster.mapping <- c("0"="Adipogenic MSC/Fib", 
                     "1"="MSC/Fib", 
                     "2"="Cycling MSC/Fib", 
                     "3"="Retinal Epi",
                     "4"="MyoFib", 
                     "5"="Early Neurons", 
                     "6"="Radial Glia", 
                     "7"="MSC/Fib",
                     "8"="Mid/Hindgut Epi", 
                     "9"="Muscle Prog", 
                     "10"="Cycling MSC/Fib", 
                     "11"="Pericytes",
                     "12"="Retinal Epi", 
                     "13"="CycProg", 
                     "14"="Retinal Epi",
                     "15"="MSC/Fib", 
                     "16"="Chondrogenic MSC/Fib", 
                     "17"="Foregut Epi", 
                     "18"="Smooth Muscle",
                     "19"="Cardiac/Skeletal Muscle", 
                     "20"="MSC/Fib", 
                     "21"="Immune", 
                     "22"="Retinal Neurons",
                     "23"="Schwann Cells", 
                     "24"="Melanoblasts", 
                     "25"="Kidney Prog", 
                     "26"="HSC")
named.ident <- plyr::revalue(ter.ident, replace = cluster.mapping)
table(named.ident)

plot.seed <- 12345
## Plot original Seurat clusters
PlotDims(umap.emb, sample.groups = ter.ident, show.axes = F, show.legend = F,
         pt.size = 0.3, alpha = 0.3, seed = plot.seed, use.brewer.pal = F,
         label.size = 4)

## Plot mapped Seurat clusters
PlotDims(umap.emb, sample.groups = named.ident, show.axes = F, show.legend = F,
         pt.size = 0.3, alpha = 0.3, seed = plot.seed, use.brewer.pal = F,
         label.size = 4)

## Save results
save.image(output.file)