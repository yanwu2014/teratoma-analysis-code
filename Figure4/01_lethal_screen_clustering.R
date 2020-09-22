library(methods)
library(Seurat)
library(swne)
library(umap)


## Set working directory to file directory
setwd("teratoma-analysis-code/Figure4")

## Parameters
output.file <- "dm-ter-screen_clustering_Seurat.Robj"
ref.seurat.file <- "../Figure1/dm-ter-human-merged_clustering_v3.seurat.Robj"

## Load the original screen dataset
screen.matrix.dir <- "../Counts/dm-ter-screen-matrices/"
screen.matrix.paths <- paste0(screen.matrix.dir, c("m1-1", "m1-2", "m2-1", "m2-2", "m3-1", "m3-2"), 
                              "/hg19/")
screen.species.paths <- paste0(screen.matrix.dir, c("m1-1", "m1-2", "m2-1", "m2-2", "m3-1", "m3-2"), 
                               "_species_class.csv")

screen.counts.list <- lapply(1:length(screen.matrix.paths), function(i) {
  counts <- Read10X(screen.matrix.paths[[i]])
  rownames(counts) <- gsub("hg19_", "", rownames(counts))
  colnames(counts) <- sapply(colnames(counts), ExtractField, field = 1, delim = "-")
    
  species.df <- read.table(screen.species.paths[[i]], sep = ",", header = T)
  human.barcodes <- as.character(subset(species.df, call == "hg19")$barcode)
  human.barcodes <- sapply(human.barcodes, ExtractField, field = 1, delim = "-")
  
  counts <- counts[,colnames(counts) %in% human.barcodes]
  colnames(counts) <- paste(colnames(counts), i, sep = "_")
  return(counts)
})
names(screen.counts.list) <- 1:length(screen.counts.list)


repool.matrix.dir <- "../Counts/dm-ter-screen-repool-matrices/"
repool.matrix.paths <- paste0(repool.matrix.dir, c("m1-1", "m1-2", "m2-1", "m2-2", "m3-1", "m3-2"))
repool.species.paths <- paste0(repool.matrix.dir, c("m1-1", "m1-2", "m2-1", "m2-2", "m3-1", "m3-2"), 
                               "_cell_class.csv")

repool.counts.list <- lapply(1:length(repool.matrix.paths), function(i) {
  counts <- Read10XHuman(repool.matrix.paths[[i]], repool.species.paths[[i]])
  colnames(counts) <- paste(colnames(counts), i + length(screen.counts.list), sep = "_")
  return(counts)
})
names(repool.counts.list) <- (length(screen.counts.list) + 1):(length(screen.counts.list) + length(repool.counts.list))


run.obj.list <- lapply(c(screen.counts.list, repool.counts.list), function(counts) {
  obj <- CreateSeuratObject(counts = counts, min.features = 200,
                            min.cells = 3)
  print(dim(obj))
  return(obj)
})

obj.list <- list()
obj.list[["ter1"]] <- merge(run.obj.list[[1]], run.obj.list[[2]])
obj.list[["ter2"]] <- merge(run.obj.list[[3]], run.obj.list[[4]])
obj.list[["ter3"]] <- merge(run.obj.list[[5]], run.obj.list[[6]])
obj.list[["ter4"]] <- merge(run.obj.list[[7]], run.obj.list[[8]])
obj.list[["ter5"]] <- merge(run.obj.list[[9]], run.obj.list[[10]])
obj.list[["ter6"]] <- merge(run.obj.list[[11]], run.obj.list[[12]])
rm(run.obj.list); gc();

obj.list <- lapply(obj.list, function(obj) {
  obj <- NormalizeData(object = obj, verbose = FALSE)
  obj <- FindVariableFeatures(object = obj, selection.method = "vst",
                              nfeatures = 3000, verbose = F)
  print(dim(obj))
  return(obj)
})

## Identify anchors and integrate (which replaces runMultiCCA)
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = 3000, dims = 1:20)
obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:20)

## Switch to integrated assay. The variable features of this assay are 
## automatically set during IntegrateData
DefaultAssay(object = obj.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
obj.integrated <- ScaleData(object = obj.integrated, features = rownames(GetAssayData(obj.integrated)))
obj.integrated <- RunPCA(object = obj.integrated, npcs = 20, verbose = F)

## Pull out 10X run
obj.integrated$run <- factor(sapply(colnames(obj.integrated), ExtractField, field = 2, delim = "_"))
names(obj.integrated$run) <- colnames(obj.integrated)
table(obj.integrated$run)

obj.integrated$teratoma <- plyr::revalue(obj.integrated$run, replace = 
                                           c("2"="1", "3"="2", "4"="2", "5"="3", "6"="3",
                                             "7"="4", "8"="4", "9"="5", "10"="5", "11"="6", "12"="6"))
levels(obj.integrated$teratoma ) <- paste0("ter", levels(obj.integrated$teratoma))
table(obj.integrated$teratoma)

## Save R object
saveRDS(obj.integrated, output.file)

## Load reference data
ref <- readRDS(ref.seurat.file)
DefaultAssay(ref) <- "integrated"

## Load reference metadata
meta.data <- read.table("../Figure1/dm-ter-human-merged_metadata.tsv", header = T, sep = "\t")
ref.ident <- meta.data$cluster; names(ref.ident) <- meta.data$cell;
table(ref.ident)

## Collapse clusters
broad.cluster.mapping <- c("Airway Epi"="Foregut Epi",
                           "Schwann Cells"="SCP",
                           "Melanoblasts"="SCP",
                           "Immune"="Hematopoietic",
                           "Erythrocyte"="Hematopoietic",
                           "HSC"="Hematopoietic",
                           "Cycling MSC/Fib"="MSC/Fib",
                           "MyoFib"="MSC/Fib",
                           "Adipogenic MSC/Fib"="MSC/Fib",
                           "Chondrogenic MSC/Fib"="MSC/Fib",
                           "Cardiac/Skeletal Muscle"="Muscle",
                           "Muscle Prog"="Muscle",
                           "Radial Glia"="Neural Prog",
                           "CycProg" = "Neural Prog",
                           "Retinal Neurons" = "Neurons",
                           "Early Neurons" = "Neurons")
ref.broad.ident <- plyr::revalue(ref.ident, replace = broad.cluster.mapping)
table(ref.broad.ident)

## Classify cells with Seurat label transfer
transfer.anchors <- FindTransferAnchors(ref, obj.integrated, dims = 1:20)
ident.df <- TransferData(transfer.anchors, refdata = ref.broad.ident, dims = 1:20)
obj.integrated$ident <- ident.df$predicted.id
obj.integrated$prediction.score.max <- ident.df$prediction.score.max

## Project data onto reference UMAP
obj.integrated <- SetIdent(obj.integrated, value = obj.integrated$ident)

ref.pc.load <- Loadings(ref, "pca")
genes.project <- intersect(rownames(ref.pc.load), rownames(GetAssayData(obj.integrated, slot = "scale.data")))
pc.emb <- t(GetAssayData(obj.integrated, slot = "scale.data")[genes.project,]) %*% ref.pc.load[genes.project,]

## Run UMAP
config <- umap.defaults
# config$metric <- "cosine"; config$min_dist <- 0.3;
umap.obj <- umap(pc.emb, config = config)
umap.emb <- umap.obj$layout

## Store in Seurat object
obj.integrated[["umap"]] <- CreateDimReducObject(embeddings = umap.emb, key = "UMAP_", assay = "integrated")

## Plot UMAP
PlotDims(umap.emb[names(obj.integrated$ident),], sample.groups = obj.integrated$ident, 
         show.axes = F, do.label = T, show.legend = F, pt.size = 0.5, alpha = 0.5, seed = 2315)

## Save R object
saveRDS(obj.integrated, output.file)
