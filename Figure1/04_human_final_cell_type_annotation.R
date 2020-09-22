library(Seurat)
library(swne)
library(cellMapper)
library(ggplot2)

## Set working directory to file directory
setwd("teratoma-analysis-code/Figure1")

## Input/output files
output.file <- "dm-ter-human-merged_mapping_v3.RData"
load(output.file)

## Update clusters with ciliated epithelium sub-clustering results
cilia.epi.clusters <- read.table("dm-ter-human-merged_cilia_epi_clusters.txt", sep = "\t",
                                 stringsAsFactors = F, header = T)
rownames(cilia.epi.clusters) <- cilia.epi.clusters$cellnames
named.ident <- as.character(named.ident); names(named.ident) <- names(ter.ident);
named.ident[rownames(cilia.epi.clusters)] <- cilia.epi.clusters$ident
named.ident <- factor(named.ident)
named.ident <- plyr::revalue(named.ident, replace = c("Airway Epithelium"="Airway Epi",
                                                      "Retinal Epithelium"="Retinal Epi"))
ter <- SetIdent(ter, value = named.ident[colnames(ter)])
table(named.ident)


## Plot original Seurat clusters
PlotDims(umap.emb, sample.groups = ter.ident, show.axes = F, show.legend = F,
         pt.size = 0.3, alpha = 0.3, seed = plot.seed, use.brewer.pal = F,
         label.size = 4)


## Set color palette and export
library(ggplot2)
plot.seed <- 3509238
g <- ggplot_build(PlotDims(umap.emb, sample.groups = named.ident, show.axes = F, show.legend = F,
                           pt.size = 0.5, alpha = 0.5, seed = plot.seed, use.brewer.pal = T,
                           do.label = F))
cluster.colors <- unique(g$data[[1]]$colour)
names(cluster.colors) <- unique(as.character(named.ident))
saveRDS(cluster.colors, file = "dm-ter-human-merged_umap_cluster_colors.Robj")


## Plot mapped clusters
pdf("dm-ter-human-merged_umap_named_clusters.pdf", width = 6, height = 6)
PlotDims(umap.emb, sample.groups = named.ident, show.axes = F, show.legend = F,
         pt.size = 0.3, alpha = 0.3, seed = plot.seed, use.brewer.pal = F,
         label.size = 4, colors.use = cluster.colors)
dev.off()

pdf("dm-ter-human-merged_umap_clusters_nolabel.pdf", width = 6, height = 6)
PlotDims(umap.emb, sample.groups = named.ident, show.axes = F, show.legend = F,
         pt.size = 0.5, alpha = 0.5, seed = plot.seed, use.brewer.pal = T,
         do.label = F, colors.use = cluster.colors)
dev.off()


## Plot by germ layer
layer.mapping <- c("MSC/Fib"="Mesoderm", 
                   "MyoFib"="Mesoderm", 
                   "Cycling MSC/Fib"="Mesoderm", 
                   "Retinal Epi"="Ectoderm", 
                   "Adipogenic MSC/Fib"="Mesoderm", 
                   "Early Neurons"="Ectoderm", 
                   "Radial Glia"="Ectoderm", 
                   "Mid/Hindgut Epi"="Endoderm", 
                   "Muscle Prog"="Mesoderm", 
                   "Pericytes"="Mesoderm",
                   "CycProg"="Ectoderm",
                   "Chondrogenic MSC/Fib"="Mesoderm", 
                   "Foregut Epi"="Endoderm",
                   "Smooth Muscle"="Mesoderm", 
                   "Cardiac/Skeletal Muscle"="Mesoderm", 
                   "Immune"="Mesoderm", 
                   "Retinal Neurons"="Ectoderm", 
                   "Schwann Cells"="Ectoderm", 
                   "Melanoblasts"="Ectoderm", 
                   "Kidney Prog"="Mesoderm", 
                   "HSC"="Mesoderm", 
                   "Airway Epi"="Endoderm",
                   "Erythrocyte"="Mesoderm")

layer.mapping.df <- data.frame(Cell_Type = names(layer.mapping), Germ_Layer = layer.mapping)
write.table(layer.mapping.df, file = "dm-ter-human-merged_cluster_layer_mapping.txt",
            sep = "\t", row.names = F, quote = F)

germ.layers <- plyr::revalue(named.ident, replace = layer.mapping)
table(germ.layers)

pdf("dm-ter-human-merged_umap_germ_layers.pdf", width = 7, height = 6)
PlotDims(umap.emb, sample.groups = germ.layers, show.axes = F, show.legend = T, do.label = F,
         pt.size = 0.5, alpha = 0.5, seed = 12535) + 
  theme(legend.title = element_blank(), legend.text = element_text(size = 14))
dev.off()

ter.id <- factor(ter$batch)
pdf("dm-ter-human-merged_umap_teratoma.pdf", width = 7, height = 6)
PlotDims(umap.emb, sample.groups = ter.id, show.axes = F, show.legend = T,
         pt.size = 0.1, alpha = 0.7, seed = 12535, do.label = F, use.brewer.pal = F)
dev.off()

barcoding <- factor(ter$barcoded)
pdf("dm-ter-human-merged_umap_barcoding.pdf", width = 7, height = 6)
PlotDims(umap.emb, sample.groups = barcoding, show.axes = F, show.legend = T,
         pt.size = 0.1, alpha = 0.7, seed = 12535, do.label = F, use.brewer.pal = F)
dev.off()


# dotplot.genes.df <- Reduce(rbind, by(ter.markers.df, ter.markers.df$cluster, head, n = 2))
# dotplot.genes.df <- Reduce(rbind, by(tf.markers.df, tf.markers.df$cluster, head, n = 2))
# genes.plot <- as.character(unique(dotplot.genes.df$gene))
key.markers.df <- read.table("dm-ter-human-merged_marker_cluster_mapping.txt", 
                             header = T, sep = "\t", stringsAsFactors = F)
genes.plot <- unique(key.markers.df$Marker)
cl.order <- unique(key.markers.df$Cell.Type)

# ter <- ScaleData(ter, genes.use = genes.plot, assay = "RNA")
heat.mat <- t(apply(GetAssayData(ter, assay = "RNA", slot = "scale.data")[genes.plot,], 1, function(x) {
  tapply(x, named.ident, mean)
}))
heat.mat[heat.mat > 5] <- 5

pdf("dm-ter-human-merged_marker_genes_heatmap.pdf", width = 6, height = 9.25)
ggHeat(heat.mat[,cl.order], clustering = "none", x.lab.size = 11, y.lab.size = 11)
dev.off()


ter$named.clusters <- factor(named.ident, levels = cl.order)
pdf("dm-ter-human-merged_marker_genes_dotplot.pdf", width = 6.5, height = 9)
DefaultAssay(ter) <- "RNA"
DotPlot(ter, features = genes.plot, group.by = "named.clusters") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_blank()) + 
  coord_flip()
dev.off()


## Export clusters and teratoma
meta.data <- data.frame(cell = colnames(ter), cluster = named.ident[colnames(ter)], 
                        teratoma = ter$batch, barcoding = ter$barcoded)
write.table(meta.data, file = "dm-ter-human-merged_metadata.tsv", sep = "\t",
            row.names = F, col.names = T, quote = F)


## Bar chart of cell type abundances in the teratoma
ident.counts <- as.integer(table(named.ident))
names(ident.counts) <- names(table(named.ident))

pdf("dm-ter-human-merged_celltypes_barplot.pdf", width = 4, height = 6.5)
ggBarplot(sort(ident.counts, decreasing = T), fill.color = "grey") + coord_flip()
dev.off()


## Save results
save.image(output.file)


## Get teratoma marker genes with named clusters
ter <- SetIdent(ter, value = named.ident[colnames(ter)])
ter.markers <- FindAllMarkers(ter, assay = "RNA", return.thresh = 0.001)
table(ter.markers$cluster)

write.table(ter.markers, file = "dm-ter-human-merged_named_cluster_markers.tsv", sep = "\t", 
            row.names = F, quote = F)


# top.ter.markers <- Reduce(rbind, by(ter.markers, ter.markers$cluster, head, n = 20))
# write.table(top.ter.markers, file = "dm-ter-human-merged_named_top_cluster_markers.tsv", sep = "\t",
#             row.names = F, quote = F)
# 
# tf.markers <- subset(ter.markers, gene %in% trrust.df$V1)
# top.tf.markers <- Reduce(rbind, by(tf.markers, tf.markers$cluster, head, n = 20))
# write.table(top.tf.markers, file = "dm-ter-human-merged_named_TF_cluster_markers.tsv", sep = "\t",
#             row.names = F, quote = F)

## Visualize marker gene expression
gene <- "PAX6"
gene.expr <- GetAssayData(ter, assay = "RNA")[gene,]
# pdf(paste0("dm-ter-human-merged_umap_", gene, "_expr.pdf"), width = 6.5, height = 6)
FeaturePlotDims(umap.emb, gene.expr, feature.name = gene, show.axes = F, pt.size = 0.5,
                alpha.plot = 0.5, quantiles = c(0.001, 0.999))
# dev.off()

## Save results
save.image(output.file)



