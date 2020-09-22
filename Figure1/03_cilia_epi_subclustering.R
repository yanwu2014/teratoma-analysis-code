library(Seurat)
library(swne)

## Set working directory to file directory
setwd("teratoma-analysis-code/Figure1")

## Save output
output.file <- "dm-ter-human-merged_cilia_epi_analysis.RData"
human.seurat.file <- "dm-ter-human-merged_clustering_v3.seurat.Robj"

## Load Seurat object
obj <- readRDS(human.seurat.file)

## Identify Seurat clusters that correspond to ciliated epithelium using marker genes and MCA mapping
## These cluster mappings may change due to randomness in the clustering algorithm
cilia.epi.clusters <- c("3", "12", "14")

## Subset to only ciliated epithelium clusters
obj <- subset(obj, ident.use = cilia.epi.clusters)
table(obj$batch)

DefaultAssay(obj) <- "RNA"
obj.list <- SplitObject(object = obj, split.by = "batch")
for(i in 1:length(obj.list)) {
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]], verbose = FALSE)
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]], selection.method = "vst", 
                                        nfeatures = 4000, verbose = F)
}

## Identify anchors and integrate (which replaces runMultiCCA)
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = 4000, k.filter = 150, 
                                      dims = 1:20)
obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:20)

## Switch to integrated assay. The variable features of this assay are 
## automatically set during IntegrateData
DefaultAssay(object = obj.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
obj.integrated <- ScaleData(obj.integrated, features = rownames(GetAssayData(obj.integrated)))
obj.integrated <- RunPCA(obj.integrated, npcs = 20, verbose = F)
ElbowPlot(obj.integrated)

obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", 
                          dims = 1:20)

## Clustering
cluster.res <- 0.1
obj.integrated <- FindNeighbors(obj.integrated, k.param = 20, reduction = "pca", dims = 1:20)
obj.integrated <- FindClusters(obj.integrated, resolution = cluster.res,
                               algorithm = 2, print.output = F)
table(Idents(obj.integrated))

## Find top markers and export for further analysis
markers.df <- FindAllMarkers(obj.integrated, assay = "RNA")
top.markers.df <- Reduce(rbind, by(subcluster.markers, subcluster.markers$cluster, head, n = 10))
write.table(top.markers.df, file = "dm-ter-human-merged_cilia_epi_subcluster_markers.tsv", sep = "\t")

## Map sub-clusters to cell types. Again these mappings may change if you're repeating the clustering
orig.mapping <- c("5"="Erythrocyte", "4"="Airway Epithelium", "0"="Retinal Epithelium",
                  "1"="Retinal Epithelium", "2"="Retinal Epithelium", "3"="Retinal Epithelium")
orig.idents <- droplevels(Idents(obj.integrated))
orig.idents <- plyr::revalue(orig.idents, replace = orig.mapping)

## Make UMAP plot
umap.emb <- Embeddings(obj.integrated, "umap")
# umap.emb <- Embeddings(obj.integrated, "tsne")

PlotDims(umap.emb, sample.groups = orig.idents, show.axes = F, show.legend = F,
         pt.size = 1.5, alpha = 0.5, seed = 43509238, use.brewer.pal = T)

batch <- obj.integrated$batch
PlotDims(umap.emb, sample.groups = batch, show.axes = F, show.legend = T,
         pt.size = 0.5, alpha = 0.5, seed = 43509238, use.brewer.pal = T,
         do.label = F)

## Save image
rm(obj, obj.list); gc();
save.image(output.file)

## Output ciliated epithelium sub-clusters
write.table(data.frame(ident = orig.idents, cellnames = names(orig.idents)), 
            file = "dm-ter-human-merged_cilia_epi_clusters.txt", row.names = F, 
            col.names = T, quote = F, sep = "\t")