#### Load datasets ####
library(Seurat)
library(swne)
library(perturbLM)
library(ggplot2)
library(Matrix)

## Set working directory to file directory
setwd("teratoma-analysis-code/Figure3")

## Output file
output.file <- "dm-ter-human-merged_brain_comparison.RData"

## Input reference files from Poliodakis et al study
fetal.cortex.matrix.1 <- "../Reference_Data/sc_dev_cortex_raw_counts_mat.RData"
fetal.cortex.metadata.1 <- "../Reference_Data/sc_dev_cortex_cell_metadata.csv"

## Load reference dataset
load(fetal.cortex.matrix.1)
metadata <- read.table(fetal.cortex.metadata.1, sep = ",", header = T)
rownames(metadata) <- metadata$Cell


obj <- CreateSeuratObject(raw_counts_mat, min.cells = 10, min.features = 200, meta.data = metadata)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- SetIdent(obj, value = "Cluster")
dim(obj)

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 4000)

ref.cluster.mapping <- c("ExDp1"="ExDp", "ExDp2"="ExDp", "ExM-U"="ExM", "InCGE"="In", "InMGE"="In",
                         "PgG2M"="CycProg", "PgS"="CycProg", "oRG"="RG", "vRG"="RG")
clusters <- plyr::revalue(Idents(obj), replace = ref.cluster.mapping)



## Load teratoma dataset and find sub-clusters
ter <- readRDS("../Figure1/dm-ter-human-merged_clustering_v3.seurat.Robj")
DefaultAssay(ter) <- "RNA"

ter.metadata <- read.table("../Figure1/dm-ter-human-merged_metadata.tsv", header = T, sep = "\t")
ter.clusters <- ter.metadata$cluster; names(ter.clusters) <- ter.metadata$cell;

ter <- SetIdent(ter, value = ter.clusters)
ter <- subset(ter, idents = c("Radial Glia", "CycProg", "Early Neurons"))

ter <- FindVariableFeatures(ter, selection.method = "vst", nfeatures = 4000)
ter <- ScaleData(ter, features = rownames(GetAssayData(ter)))
ter <- RunPCA(ter, features = VariableFeatures(ter), verbose = F)
ter <- RunUMAP(ter, dims = 1:20, verbose = F)
ter <- FindNeighbors(ter, k.param = 10)
ter <- FindClusters(ter, resolution = 0.4)

pdf("dm-ter-human-merged_neural_subclusters.pdf", width = 6, height = 5)
DimPlot(ter, reduction = "umap", label = T) + theme_void()
dev.off()


## Map teratoma sub-clusters to cell types from fetal data
var.genes <- union(VariableFeatures(ter), VariableFeatures(obj))
var.genes <- intersect(var.genes, rownames(ter))
var.genes <- intersect(var.genes, rownames(obj))

obj <- ScaleData(obj, features = var.genes)
obj <- RunPCA(obj, features = var.genes, npcs = 30, verbose = F)
obj <- FindNeighbors(obj, k.param = 10, dims = 1:20)

transfer.anchors <- FindTransferAnchors(reference = obj, query = ter, dims = 1:20, 
                                        features = var.genes)
ter.ident.df <- TransferData(anchorset = transfer.anchors, refdata = clusters, dims = 1:20)
ter$predicted.id <- ter.ident.df$predicted.id
hist(ter.ident.df$prediction.score.max)

ter.predict <- ter.ident.df
ter.predict[c("predicted.id", "prediction.score.max")] <- NULL
ter.predict <- as.matrix(ter.predict)
colnames(ter.predict) <- gsub("prediction.score.", "", colnames(ter.predict))

ter.cluster.predict <- cellMapper::SummarizeKNN(t(ter.predict), Idents(ter), assign.cells = T, 
                                                min.ident.frac = 0.1)

pdf("dm-ter-human-merged_neural_subcluster_fetal_mapping.pdf", width = 5.5, height = 3.5)
ggHeat(ter.cluster.predict, clustering = "both")
dev.off()

ter.markers <- FindAllMarkers(ter, logfc.threshold = 0.5, only.pos = T)
top.ter.markers <- do.call(rbind, by(ter.markers, ter.markers$cluster, head, n = 8))
write.table(top.ter.markers, file = "dm-ter-human-merged_neural_subcluster_markers.tsv", sep = "\t")

key.markers <- c("POU5F1", "SOX2", "HES5", "VIM", "GLI3", "HMGB2", "CCNB1", "DCX", "NEUROD1",
                 "VSX2", "PRPH", "DLX1")

DotPlot(ter, features = key.markers) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



## Map sub-clusters
ter.subcluster.mapping <- c("0"="Radial Glia", "1"="Retinal Progenitors", "2"="Early Neurons", "3"="Interneurons",
                            "4"="Peripheral Neurons", "5"="CycProg", "6"="Radial Glia", "7"=NA,
                            "8"="Early Neuroectoderm", "9"=NA)
ter.subclusters <- plyr::revalue(Idents(ter), replace = ter.subcluster.mapping)

ter$orig.cluster <- ter.clusters[colnames(ter)]
ter$sub.cluster <- factor(ter.subclusters[colnames(ter)])
ter <- SetIdent(ter, value = "sub.cluster")

DotPlot(ter, features = key.markers) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf("dm-ter-human-merged_neural_subclusters_labeled.pdf", width = 6, height = 5)
DimPlot(ter, reduction = "umap", group.by = "sub.cluster", label = T) + theme_void()
dev.off()

pdf("dm-ter-human-merged_neural_subclusters_nolabels.pdf", width = 6, height = 5)
DimPlot(ter, reduction = "umap", group.by = "sub.cluster", pt.size = 1, label = F) + theme_void() + 
  theme(legend.position = "none")
dev.off()

saveRDS(ter, file = "dm-ter-human-merged_neural_subcluster_seurat.Robj")

ter.neuro.metadata <- subset(ter.metadata, cell %in% colnames(ter))
ter.neuro.metadata$sub.cluster <- ter$sub.cluster
write.table(ter.neuro.metadata, file = "dm-ter-human-merged_neural_metadata.tsv", sep = "\t")


#### Project teratoma data onto fetal brain data ####


## Run SWNE on fetal brain data
clusters <- clusters[!is.na(clusters)]
obj <- subset(obj, cells = names(clusters))
ref.norm.counts <- ExtractNormCounts(obj)

obj <- FindNeighbors(obj, k.param = 20, dims = 1:20, prune.SNN = 0.0)

snn <- as(obj@graphs$RNA_snn, "dgCMatrix")
knn <- as(obj@graphs$RNA_nn, "dgCMatrix")
snn <- PruneSNN(snn, knn, clusters = droplevels(clusters), qval.cutoff = 1e-6)

k <- 20
n.cores <- 16 

ref.nmf.res <- RunNMF(ref.norm.counts[var.genes,], k = k, n.cores = n.cores, ica.fast = T)
ref.nmf.res$W <- ProjectFeatures(ref.norm.counts, ref.nmf.res$H, n.cores = n.cores)

ref.swne <- EmbedSWNE(ref.nmf.res$H, SNN = snn, alpha.exp = 1.25, snn.exp = 1, 
                      n_pull = 3)
ref.swne$H.coords$name <- ""

cortex.markers <- c("SOX2", ## Progenitor marker
                    "HES5", "VIM", ## RG marker
                    "HMGB2", ## Cycling progenitor marker
                    "DCX", "NEUROD1", ## Early neuron markers
                    "DLX1") ## Interneuron marker
cortex.marker.anno <- c("SOX2" = "Radial Glia", "HES5"="Radial Glia", "VIM"="Radial Glia",
                        "HMGB2"="CycProg",
                        "DCX"="Early Neurons", "NEUROD1"="Early Neurons",
                        "DLX1" = "Interneurons")
cortex.marker.anno <- factor(cortex.marker.anno)

ref.swne <- EmbedFeatures(ref.swne, ref.nmf.res$W, cortex.markers, n_pull = 3)

pdf("dm-ter-human-merged_ter_brain_swne.pdf", width = 6, height = 5)
PlotSWNE(ref.swne, alpha.plot = 0.4, sample.groups = clusters, do.label = T,
         show.legend = F, pt.size = 1, seed = 235325, use.brewer.pal = F)
dev.off()

pdf("dm-ter-human-merged_ter_brain_swne_nolabels.pdf", width = 6, height = 5)
PlotSWNE(ref.swne, alpha.plot = 0.4, sample.groups = clusters, do.label = F,
         show.legend = F, pt.size = 1, seed = 235325, use.brewer.pal = F)
dev.off()



## Project teratoma data
ter <- SetIdent(ter, value = "sub.cluster")
ter.cortex <- subset(ter, idents = c("Radial Glia", "CycProg", "Early Neurons", "Interneurons"))
ter.norm.counts <- ExtractNormCounts(ter.cortex)

ref.pc.load <- Loadings(obj, reduction = "pca")
ref.pc.emb <- t(Embeddings(obj, reduction = "pca"))

ter.cortex <- ScaleData(ter.cortex, features = var.genes)
ter.scale.counts <- GetAssayData(ter.cortex, slot = "scale.data")[var.genes,]

ter.pc.emb <- t(t(ter.scale.counts) %*% ref.pc.load[var.genes,])
ter.snn <- ProjectSNN(ter.pc.emb, ref.pc.emb, k = 10, prune.SNN = 0.0)

ter.nmf.scores <- ProjectSamples(ter.norm.counts, ref.nmf.res$W, features.use = var.genes, 
                                 n.cores = n.cores)

ter.swne <- ProjectSWNE(ref.swne, ter.nmf.scores, SNN = ter.snn, alpha.exp = 1.5,
                        snn.exp = 1, n_pull = 3)

fetal.colors <- ExtractSWNEColors(ref.swne, sample.groups = clusters, seed = 235325)
fetal.colors[["Early Neurons"]] <- fetal.colors[["ExN"]]
fetal.colors[["Interneurons"]] <- fetal.colors[["In"]]
fetal.colors[["Radial Glia"]] <- fetal.colors[["RG"]]

# pdf("dm-ter-human-merged_fetal_brain_proj_swne.pdf", width = 4, height = 4)
PlotSWNE(ter.swne, alpha.plot = 0.4, sample.groups = ter.cortex$sub.cluster, do.label = T,
         show.legend = F, pt.size = 1.75, colors.use = fetal.colors)
# dev.off()

pdf("dm-ter-human-merged_fetal_brain_proj_swne_nolabels.pdf", width = 4, height = 4)
PlotSWNE(ter.swne, alpha.plot = 0.4, sample.groups = ter.cortex$sub.cluster, do.label = F,
         show.legend = F, pt.size = 1.75, colors.use = fetal.colors)
dev.off()

FeaturePlotSWNE(ter.swne, ter.norm.counts["DLX1",])

## Save results
save.image(output.file)


## UMAP of teratoma sub-clusters using matching colors
sub.cluster.colors <- fetal.colors
sub.cluster.colors[["Early Neuroectoderm"]] <- fetal.colors[["ExM"]]
sub.cluster.colors[["Retinal Progenitors"]] <- fetal.colors[["End"]]
sub.cluster.colors[["Peripheral Neurons"]] <- fetal.colors[["Mic"]]
ter.umap.emb <- Embeddings(ter, "umap")

pdf("dm-ter-human-merged_neural_subclusters_nolabels.pdf", width = 6, height = 5)
PlotDims(ter.umap.emb, sample.groups = ter$sub.cluster, show.axes = F, alpha.plot = 0.5,
         pt.size = 1.25, show.legend = F, colors.use = sub.cluster.colors,
         do.label = F)
dev.off()


#### Correlate expression of key marker genes ####
# scale.counts <- GetAssayData(obj, slot = "scale.data")
scale.counts <- ExtractNormCounts(obj)
scale.counts <- t(apply(scale.counts[cortex.markers, names(clusters)], 1, function(x) x - mean(x)))
cluster.means <- t(apply(scale.counts, 1, function(x) tapply(x, clusters, mean)))

ter.cluster.means <- t(apply(GetAssayData(ter, slot = "scale.data")[cortex.markers,], 1, function(x) {
  tapply(x, droplevels(Idents(ter)), mean)
}))


## Radial Glia correlation
pdf("dm-ter-human-merged_ter_fetal_RG_correlation.pdf", width = 4, height = 4)
PlotCorrelation(cluster.means[cortex.markers, "RG"], ter.cluster.means[cortex.markers, "Radial Glia"],
                box = F, show.corr = T, title = "Radial Glia", x.lab = "Fetal Cortex Expression",
                y.lab = "Teratoma Expression", pt.size = 2.5, pts.label = rep(T, length(cortex.markers)),
                use.label = T, label.font.size = 5, pt.color = cortex.marker.anno) + 
  theme(legend.position = "none")
dev.off()

## INP correlation
pdf("dm-ter-human-merged_ter_fetal_INP_correlation.pdf", width = 4, height = 4)
PlotCorrelation(cluster.means[cortex.markers, "CycProg"], ter.cluster.means[cortex.markers, "CycProg"],
                box = F, show.corr = T, title = "CycProg", x.lab = "Fetal Cortex Expression",
                y.lab = "Teratoma Expression", pt.size = 2.5, pts.label = rep(T, length(cortex.markers)),
                use.label = T, label.font.size = 5, pt.color = cortex.marker.anno) + 
  theme(legend.position = "none")
dev.off()

## Early Neuron correlation
pdf("dm-ter-human-merged_ter_fetal_eNeu_correlation.pdf", width = 4, height = 4)
PlotCorrelation(cluster.means[cortex.markers, "ExN"], ter.cluster.means[cortex.markers, "Early Neurons"],
                box = F, show.corr = T, title = "Early Neurons", x.lab = "Fetal Cortex Expression",
                y.lab = "Teratoma Expression", pt.size = 2.5, pts.label = rep(T, length(cortex.markers)),
                use.label = T, label.font.size = 5, pt.color = cortex.marker.anno) + 
  theme(legend.position = "none")
dev.off()

## Interneuron correlation
pdf("dm-ter-human-merged_ter_fetal_In_correlation.pdf", width = 4, height = 4)
PlotCorrelation(cluster.means[cortex.markers, "In"], ter.cluster.means[cortex.markers, "Interneurons"],
                box = F, show.corr = T, title = "Interneurons", x.lab = "Fetal Cortex Expression",
                y.lab = "Teratoma Expression", pt.size = 2.5, pts.label = rep(T, length(cortex.markers)),
                use.label = T, label.font.size = 5, pt.color = cortex.marker.anno) + 
  theme(legend.position = "none")
dev.off()


library(cowplot)
corr.plot.legend <- get_legend(PlotCorrelation(cluster.means[cortex.markers, "ExN"], 
                                               ter.cluster.means[cortex.markers, "Early Neurons"],
                                               pt.color = cortex.marker.anno))

pdf("dm-ter-human-merged_ter_fetal_correlation_legend.pdf", width = 3, height = 3)
ggdraw(corr.plot.legend)
dev.off()



#### Plot cell type proportions ####
name.mapping <- c("ExM"="Early Neurons", "ExN"="Early Neurons", "ExDp"="Mature Neurons",
                  "In"="Interneurons", "RG"="Radial Glia")
named.clusters <- plyr::revalue(clusters, replace = name.mapping)

ter.ident.tbl <- rep(0, nlevels(named.clusters)); names(ter.ident.tbl) <- levels(named.clusters);
ter.ident.tbl[names(table(droplevels(ter.cortex$sub.cluster)))] <- table(droplevels(ter.cortex$sub.cluster))

fetal.ident.frac <- table(named.clusters)/sum(table(named.clusters))
ter.ident.frac <- (ter.ident.tbl/sum(ter.ident.tbl))[names(fetal.ident.frac)]

frac.counts <- cbind(ter.ident.frac, fetal.ident.frac)
colnames(frac.counts) <- c("Teratoma", "Fetal Brain")

frac.counts.df <- setNames(reshape2::melt(frac.counts), c("Cell_Type", "Origin", "Frac"))

pdf("dm-ter-human-merged_fetal_brain_celltypes_barplot.pdf", width = 4, height = 5)
ggplot(frac.counts.df) +
  geom_bar(aes(x = Origin, y = Frac, fill = Cell_Type), stat = "identity", 
           width = 0.75, color = "black") + 
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(hjust = 1, size = 14, angle = 90, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 15), legend.title = element_text(size = 15))
dev.off()

## Save results
save.image(output.file)




#### Developmentally stage teratoma cells ####

## Fetal cortex data from Zhong et al.
fetal.cortex.matrix.2 <- "../Reference_Data/dev_cortex_exprMatrix.tsv"
fetal.cortex.metadata.2 <- "../Reference_Data/dev_cortex_meta.tsv"

fetal.timepoint.counts <- ReadData(fetal.cortex.matrix.2)
fetal.timepoint.meta <- read.table(fetal.cortex.metadata.2, header = T, 
                                   fill = T, sep = "\t")
rownames(fetal.timepoint.meta) <- fetal.timepoint.meta$Cell
fetal.timepoint.meta <- fetal.timepoint.meta[colnames(fetal.timepoint.counts),]

obj.timepoint <- CreateSeuratObject(fetal.timepoint.counts, meta.data = fetal.timepoint.meta,
                                    min.cells = 10)
obj.timepoint <- NormalizeData(obj.timepoint)
timepoint.genes.use <- intersect(rownames(obj.timepoint), rownames(ter))

obj.timepoint <- ScaleData(obj.timepoint, features = timepoint.genes.use)
obj.timepoint <- RunPCA(obj.timepoint, features = timepoint.genes.use, npcs = 20, verbose = F)
fetal.timepoint.pc.load <- Loadings(obj.timepoint, reduction = "pca")
fetal.timepoint.pc.emb <- Embeddings(obj.timepoint, reduction = "pca")

timepoint.genes.use <- intersect(timepoint.genes.use, rownames(fetal.timepoint.pc.load))
# ter <- ScaleData(ter, features = rownames(ter))
ter.scale.data <- GetAssayData(ter, slot = "scale.data")[timepoint.genes.use,]
ter.timepoint.pc.emb <- t(ter.scale.data) %*% fetal.timepoint.pc.load[timepoint.genes.use,]

age.number <- obj.timepoint$Age_in_Weeks
age.factor <- as.character(age.number)
age.factor[age.number < 13] <- "wk5-11"
age.factor[age.number >= 13 & age.number < 14] <- "wk13"
age.factor[age.number >= 14 & age.number < 16] <- "wk14-15"
age.factor[age.number >= 16 & age.number < 18] <- "wk16-17"
age.factor[age.number >= 18] <- "wk18+"
age.factor <- factor(age.factor)

ter.expr <- rowMeans(t(ter.timepoint.pc.emb))
fetal.age.expr <- t(apply(t(fetal.timepoint.pc.emb), 1, function(x) tapply(x, age.factor, mean)))
# ter.expr <- rowMeans(ter.scale.data)
# fetal.age.expr <- t(apply(as.matrix(GetAssayData(obj.timepoint, slot = "scale.data")), 1, function(x) {
#   tapply(x, age.factor, mean)
# }))

library(gtools)
ter.age.scores <- sapply(colnames(fetal.age.expr), function(x) lsa::cosine(ter.expr, fetal.age.expr[,x]))
ter.age.scores <- ter.age.scores[mixedorder(names(ter.age.scores))]

library(perturbLM)
library(ggplot2)

pdf("dm-ter-human-merged_neural_staging_barplot.pdf", height = 3, width = 3)
ggBarplot(ter.age.scores, fill.color = "grey") + coord_flip() + 
  theme(axis.text.y = element_text(size = 16))
dev.off()

## Save results
save.image(output.file)



#### Differential expression between fetal and teratoma cell types ####

## Combine Seurat objects
Idents(ter.cortex) <- plyr::revalue(Idents(ter.cortex), replace = c("CycProg" = "Cycling Prog."))
obj <- SetIdent(obj, value = clusters[colnames(obj)])

fetal.ter.combined <- merge(obj, ter.cortex, add.cell.ids = c("fetal", "ter"))
fetal.ter.combined <- fetal.ter.combined[!grepl("MT-", rownames(fetal.ter.combined)),]
fetal.ter.combined <- NormalizeData(fetal.ter.combined)
table(Idents(fetal.ter.combined))

ter.cluster.names <- c("Radial Glia", "Cycling Prog.", "Early Neurons", "Interneurons")
fetal.cluster.names <- c("RG", "CycProg", "ExN", "In")
genes.use <- intersect(rownames(obj), rownames(ter.cortex))

## Compute logfc
ter.fetal.logfc <- lapply(1:length(ter.cluster.names), function(i) {
  fetal.expr <- rowMeans(GetAssayData(fetal.ter.combined)[, Idents(fetal.ter.combined) == fetal.cluster.names[[i]]])
  ter.expr <- rowMeans(GetAssayData(fetal.ter.combined)[, Idents(fetal.ter.combined) == ter.cluster.names[[i]]])
  log2((ter.expr + 1e-3)/(fetal.expr + 1e-3))
})
names(ter.fetal.logfc) <- ter.cluster.names


## Run GSEA
library(liger)
genesets <- ReadGenesets("go_bp_genesets.gmt")
genesets <- FilterGenesets(genesets, names(ter.fetal.logfc[[1]]), min.size = 10, max.size = 250)

ter.fetal.gsea <- MultipleGSEAEnrich(ter.fetal.logfc, genesets, n.rand = 2000, n.cores = 16)
top.ter.fetal.gsea <- subset(ter.fetal.gsea, FDR < 0.05)
top.ter.fetal.gsea <- top.ter.fetal.gsea[order(top.ter.fetal.gsea$sscore, decreasing = T),]
top.ter.fetal.gsea <- rbind(do.call(rbind, by(top.ter.fetal.gsea, top.ter.fetal.gsea$Group, head, n = 5)),
                            do.call(rbind, by(top.ter.fetal.gsea, top.ter.fetal.gsea$Group, tail, n = 5)))
top.genesets.use <- unique(top.ter.fetal.gsea$genesets)

ter.fetal.sscore <- UnflattenDataframe(ter.fetal.gsea, output.name = "sscore", row.col = "genesets",
                                       col.col = "Group")
ter.fetal.sscore <- ter.fetal.sscore[top.genesets.use,]
rownames(ter.fetal.sscore) <- gsub("GO_", "", rownames(ter.fetal.sscore))
rownames(ter.fetal.sscore) <- stringr::str_to_title(rownames(ter.fetal.sscore))
rownames(ter.fetal.sscore) <- gsub("Dna", "DNA", rownames(ter.fetal.sscore))
rownames(ter.fetal.sscore) <- gsub("rna|Rna", "RNA", rownames(ter.fetal.sscore))
rownames(ter.fetal.sscore) <- gsub("canonical_", "", rownames(ter.fetal.sscore))
rownames(ter.fetal.sscore) <- gsub("_", " ", rownames(ter.fetal.sscore))

pdf("dm-ter-human-merged_brain_ter_fetal_gsea.pdf", width = 5.25, height = 5.25)
ggHeat(ter.fetal.sscore, clustering = "row", x.lab.size = 10, y.lab.size = 10)
dev.off()


## Run differential expression
fetal.ter.de <- lapply(1:length(ter.cluster.names), function(i) {
  FindMarkers(fetal.ter.combined, ident.1 = ter.cluster.names[[i]], ident.2 = fetal.cluster.names[[i]])
})
names(fetal.ter.de) <- ter.cluster.names

top.genes.use <- unique(unlist(lapply(fetal.ter.de, function(df) {
  df <- subset(df, p_val_adj < 1e-2)
  df <- df[order(df$avg_logFC, decreasing = T),]
  c(rownames(head(df, n = 5)), rownames(tail(df, n = 5)))
})))

top.ter.fetal.logfc <- do.call(cbind, ter.fetal.logfc)[top.genes.use,]
rownames(top.ter.fetal.logfc) <- plyr::revalue(rownames(top.ter.fetal.logfc), 
                                               replace = c("ENSG00000170091"="NSG2",
                                                           "ENSG00000279576"="MALAT1"))
top.ter.fetal.logfc <- top.ter.fetal.logfc[,colnames(ter.fetal.sscore)]

pdf("dm-ter-human-merged_brain_ter_fetal_gene_logfc.pdf", width = 3, height = 5.25)
ggHeat(top.ter.fetal.logfc, clustering = "row", x.lab.size = 10, y.lab.size = 10)
dev.off()


## Save results
save.image(output.file)
