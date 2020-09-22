library(swne)
library(perturbLM)
library(Seurat)
library(Matrix)
library(glmnet)
library(ggplot2)
library(cowplot)


## Set working directory to file directory
setwd("teratoma-analysis-code/Figure4")


## Parameters
output.file <- "dm-ter-neural-screen_analysis.RData"
seurat.file <- "dm-ter-neural_clustering_Seurat.Robj"
min.cells <- 40
max.guides <- 3


#### Load clustering data and map clusters ####
ter <- readRDS(seurat.file)

# ## Optionally run de-novo clustering
# ter <- RunUMAP(ter, dims = 1:20)
# ter <- FindNeighbors(ter, dims = 1:20, k.param = 10)
# ter <- FindClusters(ter, resolution = 0.5)
# table(Idents(ter))
# 
# ## Map de-novo clusters to reference clusters
# pred.ident <- ter$predicted.id[!is.na(ter$predicted.id)]
# cl.pred.scores <- GroupOverlapCounts(UnflattenGroups(pred.ident), UnflattenGroups(Idents(ter)))
# cl.pred.scores <- apply(cl.pred.scores, 2, function(x) x/sum(x))
# ggHeat(cl.pred.scores, labRow = T)
# 
# ident.max.scores <- apply(cl.pred.scores, 2, function(x) x[which.max(x)])
# ident.mapping <- apply(cl.pred.scores, 2, function(x) names(x)[which.max(x)])
# # ident.mapping <- plyr::revalue(ident.mapping, replace = c("MSC.Fibroblast" = "MSC/Fibroblast",
# #                                                           "Cardiac.Skeletal.Muscle"="Cardiac/Skeletal Muscle",
# #                                                           "Mid.Hindgut"="Mid/Hindgut"))
# # ident.mapping <- gsub("\\.", " ", ident.mapping)
# idents <- plyr::revalue(Idents(ter), replace = ident.mapping)
# table(idents)
# 
# # idents <- idents[!is.na(idents)]
# idents.n <- as.numeric(table(idents))
# names(idents.n) <- names(table(idents))
# 
# idents <- droplevels(idents[idents %in% names(idents.n[idents.n > 200])])
# ter <- subset(ter, cells = names(idents))
# table(idents)


## Save Seurat object
idents <- factor(ter$predicted.id)
ter$ident <- idents[colnames(ter)]
# saveRDS(ter, file = seurat.file)

## Make UMAP plot
umap.emb <- Embeddings(ter, "umap")

set.seed(135235)
cell.sample <- sample(rownames(umap.emb), size = 20000)
# pdf("dm-ter-screen_proj_umap.pdf", width = 5, height = 5)
PlotDims(umap.emb[cell.sample,], sample.groups = idents[cell.sample], 
         show.axes = F, do.label = T, show.legend = F, pt.size = 0.5, 
         alpha = 0.5, seed = 31345)
# dev.off()

# # pdf("dm-ter-screen_proj_umap_nolabels.pdf", width = 5, height = 5)
# PlotDims(umap.emb[cell.sample,], sample.groups = idents[cell.sample], 
#          show.axes = F, do.label = T, show.legend = F, pt.size = 0.5, 
#          alpha = 0.5, seed = 31345, label.size = 0)
# # dev.off()

# pdf("dm-ter-screen_proj_umap.pdf", width = 5, height = 5)
PlotDims(umap.emb[cell.sample,], sample.groups = ter$run[cell.sample], 
         show.axes = F, do.label = F, show.legend = T, pt.size = 0.25, 
         alpha = 0.5, seed = 31345)
# dev.off()


#### EMD analysis ####

## Filter away NAs and small clusters
idents <- idents[!is.na(idents)]
idents <- droplevels(idents[idents %in% names(table(idents)[table(idents) > 200])])

## Load guides
min.cells <- 100
guide.dict <- ReadGroups("dm-ter-neural-screen_pheno_dict_merged.csv")
guide.dict <- lapply(guide.dict, function(x) {
  x <- gsub("\\.", "_", x)
  intersect(x, names(idents))
})
guide.dict <- guide.dict[sapply(guide.dict, length) > min.cells]
sapply(guide.dict, length)

## Remove any NTC cells that also got a true guide
ko.cells <- unique(unlist(guide.dict[!grepl("NTC", names(guide.dict))], F, F))
ntc.guide.dict <- guide.dict[grepl("NTC", names(guide.dict))]
ntc.guide.dict <- lapply(ntc.guide.dict, function(x) x[!x %in% ko.cells])
guide.dict[names(ntc.guide.dict)] <- ntc.guide.dict
guide.dict[["NTC"]] <- unique(unlist(ntc.guide.dict, F, F))
sapply(guide.dict, length)


## Make gene dictionary
genes <- unique(sapply(names(guide.dict), ExtractField, field = 1, delim = "-"))
gene.dict <- lapply(genes, function(g) {
  unique(unlist(guide.dict[grepl(g, names(guide.dict))]))
})
names(gene.dict) <- genes
sapply(gene.dict, length)


## Compute cluster distances
umap.emb <- Embeddings(ter, "umap")[names(idents),]
cl.emb <- apply(t(umap.emb), 1, function(x) tapply(x, idents, mean))
cl.dist <- as.matrix(dist(cl.emb))

## Compute guide level EMD
guide.emd <- computeEMD(guide.dict, idents, cl.dist)
guide.emd.ntc <- sort(guide.emd["NTC",])
ggBarplot(guide.emd.ntc, fill.color = "skyblue")

## Compute gene level EMD
gene.emd <- computeEMD(gene.dict, idents, cl.dist)
gene.emd.ntc <- sort(gene.emd["NTC",])
ggBarplot(gene.emd.ntc, fill.color = "skyblue")

## Compare gene level EMD with embryonic lethal screen EMD
screen.avg.emd <- readRDS("dm-ter-screen_avg_emd.Robj")
screen.avg.emd <- screen.avg.emd[c("NTC1", "NTC2", "NTC3", "NTC4", "NTC5",
                                   "KLF6", "CDX2", "RUNX1", "ASCL1", "TWIST1")]
gene.emd.ntc <- gene.emd.ntc[names(gene.emd.ntc) != "NTC"]

## Normalize NTC EMDs to have an average of 1
neural.scale.factor <- 1/mean(gene.emd.ntc[grepl("NTC", names(gene.emd.ntc))])
screen.scale.factor <- 1/mean(screen.avg.emd[grepl("NTC", names(screen.avg.emd))])

pdf("dm-ter-neural_EMD_comparison_barplot_normalized.pdf", width = 6, height = 3)
ggBarplot(c(screen.avg.emd*screen.scale.factor, gene.emd.ntc*neural.scale.factor), 
          fill.color = "skyblue")
dev.off()


## Save results
save.image(output.file)

