library(methods)
library(Seurat)
library(swne)
library(ggplot2)

## Set working directory to file directory
setwd("teratoma-analysis-code/Figure2")

## Output file
output.file <- "dm-ter-chimera-merged_analysis.RData"

## Paths to merged cell line counts matrix, species identification, doublet detection files
## These paths may change depending on where you download your data to
input.counts.file <- "../Counts/dm-ter-chimera-merged/"
input.species.file <- "../Counts/dm-ter-chimera-merged_species_class.csv"
input.doublets.file <- "../Counts/dm-ter-chimera-merged_doublets.txt"

## Load the dataset
counts <- Read10X(input.counts.file)
counts <- counts[grepl("hg19", rownames(counts)),]

cell.class.df <- read.table(input.species.file, sep = ",", header = T,
                            stringsAsFactors = F)
cell.class.df <- subset(cell.class.df, call == "hg19")
counts <- counts[,colnames(counts) %in% cell.class.df$barcode]

counts <- counts[Matrix::rowSums(counts) > 20,]; dim(counts);
# write.table(as.matrix(counts), file = "dm-ter-chimera-merged.dense.matrix.tsv", 
#             sep = "\t")

## Remove doublets in barcoded dataset filtered with DoubletDetection
doublets <- as.logical(scan(input.doublets.file, sep = "\n", what = numeric()))
stopifnot(length(doublets) == ncol(counts))
doublets[is.na(doublets)] <- 1; sum(doublets)/length(doublets);
counts <- counts[,!doublets]
dim(counts)

## Remove noncoding RNAs and mitochondrial genes
rownames(counts) <- sapply(rownames(counts), function(x) gsub("hg19_", "", x))

## Pull out teratoma id
batch <- factor(sapply(colnames(counts), ExtractField, field = 2, delim = "\\-"))
names(batch) <- colnames(counts)
levels(batch) <- c("H9", "HUES62", "PGP1")

## Create and setup seurat objects for each teratoma
obj <- CreateSeuratObject(counts = counts, meta.data = data.frame(batch), min.features = 200,
                          min.cells = round(0.001*ncol(counts)))
obj.list <- SplitObject(object = obj, split.by = "batch"); rm(obj); invisible(gc());
for(i in 1:length(obj.list)) {
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]], verbose = FALSE)
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]], selection.method = "vst", 
                                        nfeatures = 4000, verbose = F)
}

## Identify anchors and integrate (which replaces runMultiCCA)
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = 3000, dims = 1:30)
obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:30)

## Switch to integrated assay. The variable features of this assay are 
## automatically set during IntegrateData
DefaultAssay(object = obj.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
obj.integrated <- ScaleData(object = obj.integrated, features = rownames(GetAssayData(obj.integrated)))
obj.integrated <- RunPCA(object = obj.integrated, features = rownames(GetAssayData(obj.integrated)),
                         verbose = F)
ElbowPlot(obj.integrated, ndims = 50)

obj.integrated <- RunUMAP(object = obj.integrated, reduction = "pca", 
                          dims = 1:30)


## Map cells to H1 teratoma clusters
ref.seurat.file <- "../Figure1/dm-ter-human-merged_clustering_v3.seurat.Robj"
ref <- readRDS(ref.seurat.file)
DefaultAssay(ref) <- "integrated"

## Load metadata
meta.data <- read.table("../Figure1/dm-ter-human-merged_metadata.tsv", header = T, sep = "\t")
ref.ident <- meta.data$cluster; names(ref.ident) <- meta.data$cell;
table(ref.ident)

## Classify cells with Seurat label transfer
transfer.anchors <- FindTransferAnchors(reference = ref, obj.integrated, dims = 1:30)
chimera.ident.df <- TransferData(anchorset = transfer.anchors, refdata = ref.ident, dims = 1:30)
chimera.ident <- chimera.ident.df$predicted.id; names(chimera.ident) <- rownames(chimera.ident.df);
table(chimera.ident)

## Write cell line metadata to file
chimera.meta.df <- data.frame(cell = names(chimera.ident),
                              ident = chimera.ident,
                              cell_line = batch[names(chimera.ident)])
write.table(chimera.meta.df, file = "dm-ter-chimera-merged_metadata.tsv",
            sep = "\t", quote = F)

## Project new cell lines onto H1 umap
library(umap)

ref.pc.load <- Loadings(ref, "pca")
genes.project <- intersect(rownames(ref.pc.load), rownames(GetAssayData(obj.integrated, slot = "scale.data")))
pc.emb <- t(GetAssayData(obj.integrated, slot = "scale.data")[genes.project,]) %*% ref.pc.load[genes.project,]

umap.cfg <- umap.defaults
umap.cfg$min_dist <- 0.3; umap.cfg$metric <- "cosine";
umap.obj <- umap(pc.emb, config = umap.cfg)
umap.emb <- umap.obj$layout

## Named cluster plots
plot.seed <- 3509238
cluster.colors <- readRDS("../Figure1/dm-ter-human-merged_umap_cluster_colors.Robj")

pdf("dm-ter-chimera-merged_named_cluster_umap.pdf", width = 8, height = 6)
PlotDims(umap.emb, sample.groups = chimera.ident, show.legend = T, show.axes = F, do.label = T,
         label.size = 4, alpha = 0.5, pt.size = 0.3, colors.use = cluster.colors)
dev.off()

pdf("dm-ter-chimera-merged_named_cluster_umap_nolabels.pdf", width = 6, height = 6)
PlotDims(umap.emb, sample.groups = chimera.ident, show.legend = F, show.axes = F, do.label = F,
         alpha = 0.5, pt.size = 0.3, colors.use = cluster.colors)
dev.off()

pdf("dm-ter-chimera-merged_teratoma_umap.pdf", width = 7, height = 6)
PlotDims(umap.emb, sample.groups = batch, show.legend = T, show.axes = F, do.label = F,
         label.size = 4, alpha = 0.3, pt.size = 0.2, seed = plot.seed, use.brewer.pal = F)
dev.off()

## Save results
save.image(output.file)

