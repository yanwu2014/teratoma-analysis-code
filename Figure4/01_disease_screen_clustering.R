library(methods)
library(Seurat)
library(swne)
library(Matrix)
library(umap)


## Set working directory to file directory
setwd("teratoma-analysis-code/Figure4")

## Parameters
output.file <- "dm-ter-neural_clustering_Seurat.Robj"
ref.seurat.file <- "../Figure1/dm-ter-human-merged_clustering_v3.seurat.Robj"


## Load the original screen dataset
matrix.dir <- "../Counts/dm-ter-neural-matrices/"
samples <- c("dm-ter-neural-m1-1", 
             "dm-ter-neural-m1-2", 
             "dm-ter-neural-m2-1", 
             "dm-ter-neural-m2-2")
matrix.paths <- paste0(matrix.dir, samples)
species.paths <- paste0(matrix.dir, samples, "_species_class.csv")

counts.list <- lapply(1:length(matrix.paths), function(i) {
  counts <- Read10X(matrix.paths[[i]])
  rownames(counts) <- gsub("hg19_", "", rownames(counts))
  counts <- counts[!grepl("mm10", rownames(counts)),]
  
  species.df <- read.table(species.paths[[i]], sep = ",", header = T)
  human.barcodes <- as.character(subset(species.df, call == "hg19")$barcode)
  # human.barcodes <- sapply(human.barcodes, ExtractField, field = 1, delim = "-")
  
  counts <- counts[,colnames(counts) %in% human.barcodes]
  colnames(counts) <- sapply(colnames(counts), ExtractField, field = 1, delim = "-")
  colnames(counts) <- paste(colnames(counts), i, sep = "_")
  return(counts)
})
names(counts.list) <- 1:length(counts.list)


run.obj.list <- lapply(counts.list, function(counts) {
  obj <- CreateSeuratObject(counts = counts, min.features = 200,
                            min.cells = 3)
  print(dim(obj))
  return(obj)
})

for(i in 1:length(run.obj.list)) {
  out.bc.file <- paste0(matrix.paths[[i]], ".human.barcodes.tsv")
  out.bc <- sapply(colnames(run.obj.list[[i]]), ExtractField, field = 1, delim = "_")
  write(out.bc, file = out.bc.file, sep = "\n")
}

obj.list <- list()
obj.list[["ter1"]] <- merge(run.obj.list[[1]], run.obj.list[[2]])
obj.list[["ter2"]] <- merge(run.obj.list[[3]], run.obj.list[[4]])
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

## Scale data and run PCA
obj.integrated <- ScaleData(object = obj.integrated, features = rownames(GetAssayData(obj.integrated)))
obj.integrated <- RunPCA(object = obj.integrated, npcs = 30, verbose = F)
ElbowPlot(obj.integrated, ndims = 30)

## Run UMAP
obj.integrated <- RunUMAP(object = obj.integrated, dims = 1:20, reduction = "pca", verbose = F)

## Pull out 10X run
obj.integrated$run <- factor(sapply(colnames(obj.integrated), ExtractField, field = 2, delim = "_"))
names(obj.integrated$run) <- colnames(obj.integrated)
table(obj.integrated$run)

obj.integrated$teratoma <- plyr::revalue(obj.integrated$run, replace = c("2"="1", "3"="2", "4"="2"))
levels(obj.integrated$teratoma ) <- paste0("ter", levels(obj.integrated$teratoma))
table(obj.integrated$teratoma)

umap.emb <- Embeddings(obj.integrated, "umap")
PlotDims(umap.emb, obj.integrated$teratoma, show.axes = F, show.legend = T, do.label = F,
         pt.size = 0.1)


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
ref.tissue.ident <- plyr::revalue(ref.ident, replace = broad.cluster.mapping)
table(ref.tissue.ident)

## Classify cells with Seurat label transfer
transfer.anchors <- FindTransferAnchors(ref, obj.integrated, dims = 1:20)
ident.df <- TransferData(transfer.anchors, refdata = ref.tissue.ident, dims = 1:20)
ident.df$predicted.id[ident.df$prediction.score.max < 0.5] <- NA

# obj.integrated@meta.data[grepl("prediction.score", colnames(obj.integrated@meta.data))] <- NULL
obj.integrated <- AddMetaData(obj.integrated, metadata = ident.df)
hist(ident.df$prediction.score.max)

# ## Project data onto reference UMAP
# obj.integrated <- SetIdent(obj.integrated, value = obj.integrated$ident)
# 
# ref.pc.load <- Loadings(ref, "pca")
# genes.project <- intersect(rownames(ref.pc.load), rownames(GetAssayData(obj.integrated, slot = "scale.data")))
# pc.emb <- t(GetAssayData(obj.integrated, slot = "scale.data")[genes.project,]) %*% ref.pc.load[genes.project,]
# 
# ## Run UMAP
# config <- umap.defaults
# # config$metric <- "cosine"; config$min_dist <- 0.3;
# umap.obj <- umap(pc.emb, config = config)
# umap.emb <- umap.obj$layout
# 
# ## Store in Seurat object
# obj.integrated[["umap"]] <- CreateDimReducObject(embeddings = umap.emb, key = "UMAP_", assay = "integrated")
# 
# ## Plot UMAP
# PlotDims(umap.emb[names(obj.integrated$ident),], sample.groups = obj.integrated$ident, 
#          show.axes = F, do.label = T, show.legend = F, pt.size = 0.5, alpha = 0.5, seed = 2315)
# PlotDims(umap.emb[names(obj.integrated$ident),], sample.groups = obj.integrated$teratoma, 
#          show.axes = F, do.label = F, show.legend = T, pt.size = 0.1, alpha = 0.5, seed = 2315)

## Save R object
# saveRDS(obj.integrated, output.file)
saveRDS(obj.integrated, "dm-ter-neural_clustering_Seurat_noCollapse.Robj")
