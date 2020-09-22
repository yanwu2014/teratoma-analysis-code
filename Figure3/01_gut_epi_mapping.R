## Load reference data
library(Seurat)
library(swne)
library(perturbLM)
library(ggplot2)
library(gtools)
library(Matrix)

## Set working directory to file directory
setwd("teratoma-analysis-code/Figure3")

## Output file to save to
output.file <- "dm-ter-human-merged_gut_comparison.RData"


#### Load data ####
fetal.counts <- ReadData("../Reference_Data/fetal_gut_counts.txt")
fetal.metadata <- read.table("../Reference_Data/fetal_gut_metadata.txt", sep = "\t", header = T)
fetal.metadata$Sample <- sapply(as.character(fetal.metadata$Sample), function(x) gsub("W", "w", x))
rownames(fetal.metadata) <- fetal.metadata$Sample

fetal.obj <- CreateSeuratObject(fetal.counts, min.cells = 10, min.features = 200, 
                              meta.data = fetal.metadata)
fetal.obj <- SetIdent(fetal.obj, value = "CellType")
fetal.obj <- subset(fetal.obj, idents = "Epithelial")

age <- sapply(colnames(fetal.obj), ExtractField, field = 2, delim = "_")
fetal.obj <- AddMetaData(fetal.obj, age, col.name = "Age")
table(age)

# fetal.obj <- SubsetData(fetal.obj, ident.use = c("6w", "7w", "8w", "9w", "11w"))

fetal.obj <- NormalizeData(fetal.obj)
fetal.obj <- FindVariableFeatures(fetal.obj, nfeatures = 1e3)

fetal.obj <- ScaleData(fetal.obj, features = rownames(fetal.obj))
fetal.obj <- RunPCA(fetal.obj, features = VariableFeatures(fetal.obj), verbose = F)
fetal.obj <- FindNeighbors(fetal.obj, dims = 1:20, k.param = 30, prune.SNN = 0)
fetal.obj <- RunUMAP(fetal.obj, dims = 1:20, verbose = F)

fetal.norm.counts <- ExtractNormCounts(fetal.obj)
fetal.var.genes <- VariableFeatures(fetal.obj)
# fetal.var.genes <- rownames(fetal.obj)

## Get cell types
fetal.ident.mapping <- c("LI"="Mid/Hindgut Epi", "SI"="Mid/Hindgut Epi", "stomach"="Foregut Epi",
                         "esophagus"="Foregut Epi")
fetal.ident <- plyr::revalue(fetal.obj$Tissue, replace = fetal.ident.mapping)



#### Run SWNE on fetal data ####
k <- 12
n.cores <- 16 

fetal.nmf.res <- RunNMF(fetal.norm.counts[fetal.var.genes,], k = k, n.cores = n.cores, ica.fast = T)
fetal.nmf.res$W <- ProjectFeatures(fetal.norm.counts, fetal.nmf.res$H, n.cores = n.cores)

knn <- as(fetal.obj@graphs$RNA_nn, "dgCMatrix")
snn <- as(fetal.obj@graphs$RNA_snn, "dgCMatrix")
snn <- PruneSNN(snn, knn, clusters = , qval.cutoff = 1e-3)

fetal.swne <- EmbedSWNE(fetal.nmf.res$H, SNN = snn, alpha.exp = 1, snn.exp = 0.25, n_pull = 3)
fetal.swne$H.coords$name <- ""

features.embed <- c("PAX9", "SOX2", "CDX2", "HHEX", "CDX1", "FOXJ1")
fetal.swne <- EmbedFeatures(fetal.swne, fetal.nmf.res$W, features.embed, n_pull = 3)

pdf("dm-ter-human-merged_ter_fetal_swne.pdf", width = 4, height = 4)
PlotSWNE(fetal.swne, alpha.plot = 0.4, sample.groups = fetal.ident, do.label = F,
         show.legend = F, pt.size = 1.5, seed = 1235325)
dev.off()



## Load teratoma data
obj <- readRDS("../Figure1/dm-ter-human-merged_clustering_v3.seurat.Robj")
ter.metadata <- read.table("../Figure1/dm-ter-human-merged_metadata.tsv", header = T, sep = "\t")
ter.clusters <- ter.metadata$cluster; names(ter.clusters) <- ter.metadata$cell;

## Subset to gut epithelium
obj <- SetIdent(obj, value = ter.clusters[colnames(obj)])
obj <- subset(obj, idents = c("Foregut Epi", "Mid/Hindgut Epi"))

DefaultAssay(obj) <- "RNA"
obj <- ScaleData(obj, features = rownames(obj))

## Project teratoma data onto fetal data
genes.use <- intersect(VariableFeatures(fetal.obj), rownames(obj))

fetal.pc.emb <- Embeddings(fetal.obj, "pca")[,1:20]
fetal.pc.load <- Loadings(fetal.obj, "pca")[,1:20]
ter.pc.emb <- t(t(GetAssayData(obj, slot = "scale.data")[genes.use,]) %*% fetal.pc.load[genes.use,])

ter.ident <- droplevels(plyr::revalue(Idents(obj), replace = c("8"="Mid/Hindgut Epi", "17"="Foregut Epi")))
table(ter.ident)

ter.norm.counts <- ExtractNormCounts(obj)
ter.nmf.scores <- ProjectSamples(ter.norm.counts, fetal.nmf.res$W, features.use = genes.use, n.cores = n.cores)
ter.snn <- ProjectSNN(ter.pc.emb, t(fetal.pc.emb), k = 20, prune.SNN = 0)

## Project fetal data onto teratoma data
ter.swne <- ProjectSWNE(fetal.swne, ter.nmf.scores, SNN = ter.snn, alpha.exp = 1,
                        snn.exp = 1, n_pull = 3)

colors.use <- ExtractSWNEColors(fetal.swne, fetal.ident, 1235325)

pdf("dm-ter-human-merged_ter_gut_proj_swne.pdf", width = 4, height = 4)
PlotSWNE(ter.swne, alpha.plot = 0.4, sample.groups = ter.ident, do.label = F,
         show.legend = F, pt.size = 1.5, colors.use = colors.use)
dev.off()

gene <- "SPRR3"
FeaturePlotSWNE(ter.swne, ter.norm.counts[gene,], feature.name = gene, alpha.plot = 0.4,
                pt.size = 1.5, quantiles = c(0.005, 0.995))



#### Correlate key marker genes between teratoma and fetal data ####
ter.cluster.means <- t(apply(GetAssayData(obj, slot = "scale.data"), 1, function(x) {
  tapply(x, ter.ident, mean)
}))
cluster.means <- t(apply(GetAssayData(fetal.obj, slot = "scale.data"), 1, function(x) {
  tapply(x, fetal.ident, mean)
}))

genes.corr <- c("PAX9", "SOX2", "CDX2", "CDX1", "HHEX", "FOXJ1")
genes.corr.mapping <- c("PAX9"="Foregut", "SOX2"="Foregut", "FOXJ1"="Foregut Epi", 
                        "HHEX"="Mid/Hindgut Epi",
                        "CDX1"="Mid/Hindgut Epi", "CDX2"="Mid/Hindgut Epi")

## Mid/Hindgut correlation
pdf("dm-ter-human-merged_ter_fetal_mid_hindgut_correlation.pdf", width = 4, height = 4)
PlotCorrelation(cluster.means[genes.corr, "Mid/Hindgut Epi"], ter.cluster.means[genes.corr, "Mid/Hindgut Epi"],
                box = F, show.corr = T, title = "Mid/Hindgut Epi", x.lab = "Fetal Gut Expression",
                y.lab = "Teratoma Expression", pt.size = 2.5, pts.label = rep(T, length(genes.corr)),
                use.label = T, label.font.size = 5, pt.color = genes.corr.mapping) + 
  theme(legend.position = "none")
dev.off()

## Foregut correlation
pdf("dm-ter-human-merged_ter_fetal_foregut_correlation.pdf", width = 4, height = 4)
PlotCorrelation(cluster.means[genes.corr, "Foregut Epi"], ter.cluster.means[genes.corr, "Foregut Epi"],
                box = F, show.corr = T, title = "Foregut Epi", x.lab = "Fetal Gut Expression",
                y.lab = "Teratoma Expression", pt.size = 2.5, pts.label = rep(T, length(genes.corr)),
                use.label = T, label.font.size = 5, pt.color = genes.corr.mapping) + 
  theme(legend.position = "none")
dev.off()

library(cowplot)
corr.plot.legend <- get_legend(PlotCorrelation(cluster.means[genes.corr, "Foregut Epi"], 
                                               ter.cluster.means[genes.corr, "Foregut Epi"],
                                               pt.color = genes.corr.mapping))

pdf("dm-ter-human-merged_ter_fetal_gut_correlation_legend.pdf", width = 3, height = 3)
ggdraw(corr.plot.legend)
dev.off()



#### Stage teratoma data ####
ter.expr <- rowMeans(ter.pc.emb)
fetal.age.expr <- t(apply(t(fetal.pc.emb), 1, function(x) tapply(x, fetal.obj$Age, mean)))

ter.age.scores <- sapply(colnames(fetal.age.expr), function(x) lsa::cosine(ter.expr, fetal.age.expr[,x]))
ter.age.scores <- ter.age.scores[mixedorder(names(ter.age.scores))]

pdf("dm-ter-human-merged_gut_staging_barplot.pdf", height = 4, width = 3)
ggBarplot(ter.age.scores, fill.color = "grey") + coord_flip() + 
  theme(axis.text.y = element_text(size = 16))
dev.off()



#### Plot cell type proportions ####
ter.ident.frac <- table(ter.ident)/sum(table(ter.ident))
fetal.ident.frac <- table(fetal.ident)/sum(table(fetal.ident))
frac.counts <- cbind(ter.ident.frac, fetal.ident.frac)
colnames(frac.counts) <- c("Teratoma", "Fetal Gut")

frac.counts.df <- setNames(reshape2::melt(frac.counts), c("Cell_Type", "Origin", "Frac"))

pdf("dm-ter-human-merged_fetal_gut_celltypes_barplot.pdf", width = 5.5, height = 3)
ggplot(frac.counts.df) +
  geom_bar(aes(x = Origin, y = Frac, fill = Cell_Type), stat = "identity", 
           width = 0.75, color = "black") + 
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(hjust = 1, size = 14, angle = 90, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16)) + 
  scale_fill_brewer(palette = "Set2") + coord_flip()
dev.off()


## Save results
save.image(output.file)


#### Differential expression between fetal and teratoma cell types ####

## Combine Seurat objects
obj <- SetIdent(obj, value = ter.ident[colnames(obj)])
levels(Idents(obj)) <- c("ter.Mid/Hindgut Epi", "ter.Foregut Epi")
fetal.obj <- SetIdent(fetal.obj, value = fetal.ident[colnames(fetal.obj)])

fetal.ter.combined <- merge(fetal.obj, obj, add.cell.ids = c("fetal", "ter"))
fetal.ter.combined <- fetal.ter.combined[!grepl("MT-", rownames(fetal.ter.combined)),]
fetal.ter.combined <- NormalizeData(fetal.ter.combined)
table(Idents(fetal.ter.combined))

ter.cluster.names <- c("ter.Foregut Epi", "ter.Mid/Hindgut Epi")
fetal.cluster.names <- c("Foregut Epi", "Mid/Hindgut Epi")

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
top.ter.fetal.gsea <- subset(ter.fetal.gsea, FDR < 0.1)
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
rownames(ter.fetal.sscore) <- gsub("canonical_|containing_", "", rownames(ter.fetal.sscore))
rownames(ter.fetal.sscore) <- gsub("_", " ", rownames(ter.fetal.sscore))
colnames(ter.fetal.sscore) <- gsub("ter.", "", colnames(ter.fetal.sscore))

pdf("dm-ter-human-merged_gut_ter_fetal_gsea.pdf", width = 3.5, height = 4)
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
colnames(top.ter.fetal.logfc) <- gsub("ter.", "", colnames(top.ter.fetal.logfc))
top.ter.fetal.logfc <- top.ter.fetal.logfc[,colnames(ter.fetal.sscore)]

pdf("dm-ter-human-merged_gut_ter_fetal_gene_logfc.pdf", width = 2.25, height = 4)
ggHeat(top.ter.fetal.logfc, clustering = "row", x.lab.size = 10, y.lab.size = 10)
dev.off()


## Save results
save.image(output.file)
