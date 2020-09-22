library(methods)
library(Seurat)
library(swne)

## Set working directory to file directory
setwd("teratoma-analysis-code/Figure1")

## Paths to input/output files. 
## These may change depending on where you download your data to
output.file <- "dm-ter-mouse-merged_clustering_v3.RData"
input.counts.file <- "../Counts/dm-ter-mouse-merged/"
input.doublets.file <- "../Counts/dm-ter-mouse-merged_doublets.txt"

## Load the dataset
counts <- ReadData()

## Remove doublets in barcoded dataset filtered with DoubletDetection
doublets <- as.logical(scan(input.doublets.file, sep = "\n", what = numeric()))
doublets[is.na(doublets)] <- 1; sum(doublets)/length(doublets);
counts <- counts[,!doublets]

## Remove noncoding RNAs and mitochondrial genes
rownames(counts) <- sapply(rownames(counts), function(x) gsub("hg19_", "", x))
rownames(counts) <- sapply(rownames(counts), function(x) gsub("mm10_", "", x))
genes.keep <- rownames(counts)[!grepl("RP11|MT", rownames(counts))]
counts <- counts[genes.keep,]

## Pull out teratoma id
batch <- factor(sapply(colnames(counts), ExtractField, field = 2, delim = "\\."))
names(batch) <- colnames(counts)

## Pull out barcoded vs non-barcoded
barcoded <- as.character(batch)
barcoded[barcoded %in% c("1", "2", "3", "4")] <- "nonbc"
barcoded[barcoded %in% c("5", "6", "7")] <- "bc"
barcoded <- factor(barcoded); names(barcoded) <- colnames(counts);

## Rename batch IDs
levels(batch) <- paste0("ter", levels(batch))

## Create and setup seurat objects for each teratoma
meta.data <- data.frame(batch, barcoded)
obj <- CreateSeuratObject(counts = counts, meta.data = meta.data, min.features = 200,
                          min.cells = round(0.001*ncol(counts)))
obj.list <- SplitObject(object = obj, split.by = "batch"); rm(obj); invisible(gc());
for(i in 1:length(obj.list)) {
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]], verbose = FALSE)
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]], selection.method = "vst", 
                                        nfeatures = 4000, verbose = F)
}

## Identify anchors and integrate (which replaces runMultiCCA)
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = 4000, dims = 1:30)
obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:30)

## Switch to integrated assay. The variable features of this assay are 
## automatically set during IntegrateData
DefaultAssay(object = obj.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
obj.integrated <- ScaleData(object = obj.integrated, features = rownames(GetAssayData(obj.integrated)))
obj.integrated <- RunPCA(object = obj.integrated, npcs = 30, verbose = F)
obj.integrated <- RunUMAP(object = obj.integrated, reduction = "pca", 
                          dims = 1:30, metric = "cosine")

## Clustering
cluster.res <- 0.6
obj.integrated <- FindNeighbors(obj.integrated, k.param = 30, reduction = "pca", dims = 1:30)
obj.integrated <- FindClusters(obj.integrated, resolution = cluster.res, 
                               algorithm = 2, print.output = F)
table(Idents(obj.integrated))

DimPlot(obj.integrated, reduction = "umap", label = T)

markers.df <- FindAllMarkers(obj.integrated, assay = "RNA", only.pos = F, print.bar = T)
table(markers.df$cluster)

## Save image
save.image(output.file)

## Save Seurat object
saveRDS(obj.integrated, file = gsub(".RData", ".seurat.Robj", output.file))

## Save marker genes
write.table(markers.df, file = gsub(".RData", ".markers.tsv", output.file), sep = "\t")
