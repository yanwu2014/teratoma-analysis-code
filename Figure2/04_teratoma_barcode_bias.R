## Cell type variability analysis
library(Seurat)
library(swne)
library(perturbLM)
library(ggplot2)
library(entropy)
library(cowplot)


## Set working directory to file directory
setwd("teratoma-analysis-code/Figure2")

## Output file
output.file <- "dm-ter-human-merged_teratoma_barcoding_analysis.RData"


## Load cell type and teratoma assignments
meta.data <- read.table("../Figure1/dm-ter-human-merged_metadata.tsv", header = T, 
                        row.names = NULL, sep = "\t")

clusters <- meta.data$cluster; names(clusters) <- meta.data$cell;
clusters.list <- UnflattenGroups(clusters)

ter.id <- meta.data$teratoma; names(ter.id) <- meta.data$cell;
ter.id.list <- UnflattenGroups(ter.id)

layer.mapping.df <- read.table("../Figure1/dm-ter-human-merged_cluster_layer_mapping.txt", header = T, sep = "\t")
layer.mapping <- layer.mapping.df$Germ_Layer; names(layer.mapping) <- layer.mapping.df$Cell_Type;

ter.cluster.counts <- GroupOverlapCounts(clusters.list, ter.id.list)
cluster.div <- apply(ter.cluster.counts, 1, function(x) {
  KL.empirical(x, colSums(ter.cluster.counts))*sum(x)
})


#### Barcode enrichment analysis ####
barcode.list <- ReadGroups("dm-ter-bc-human_pheno_dict.csv")
barcode.list <- lapply(barcode.list, function(x) {
  tag <- sapply(x, ExtractField, field = 1, delim = "\\.")
  id <- as.numeric(sapply(x, ExtractField, field = 2, delim = "\\.")) + 4
  paste(tag, id, sep = ".")
})

## Look at fraction of barcoded cells in each teratoma/cluster
cells <- as.character(subset(meta.data, teratoma %in% c("ter5", "ter6", "ter7"))$cell)
ter <- droplevels(subset(meta.data, teratoma %in% c("ter5", "ter6", "ter7"))$teratoma)
names(ter) <- cells
ter.list <- UnflattenGroups(ter)

barcoded.cells <- unique(unlist(barcode.list, F, F))
ter.count <- sapply(ter.list, length)
ter.barcode.count <- sapply(ter.list, function(x) sum(x %in% barcoded.cells))
ter.barcode.frac <- ter.barcode.count/ter.count
ter.sc.barcode.df <- data.frame(n_cells = ter.count, bc_cells = ter.barcode.count, frac = ter.barcode.frac)
write.table(ter.sc.barcode.df, file = "dm-ter-bc-human_single_cell_stats.txt", sep = "\t")


## Compute whether certain clusters are seen more often in certain barcodes and vice versa
barcode.cluster.counts.raw <- GroupOverlapCounts(barcode.list, clusters.list)
cluster.unique.barcodes <- apply(barcode.cluster.counts.raw, 2, function(x) sum(x > 0))
cluster.barcode.cells <- colSums(barcode.cluster.counts.raw)
cluster.barcode.frac <- sort(cluster.unique.barcodes/cluster.barcode.cells)



## Require at least 10 cells per barcode when computing KL divergence for barcodes
barcode.cluster.counts <- barcode.cluster.counts.raw[rowSums(barcode.cluster.counts.raw) > 10,]
dim(barcode.cluster.counts)

barcode.div <- apply(barcode.cluster.counts, 1, function(x) {
  KL.empirical(x, colSums(barcode.cluster.counts))*sum(x)
})
barcode.div <- sort(barcode.div, decreasing = T)

library(ggplot2)
pdf("dm-ter-bc_barcode_bias_barplot.pdf", width = 6, height = 2)
ggBarplot(barcode.div, fill.color = "lightgrey") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
dev.off()

cluster.barcode.div <- apply(barcode.cluster.counts, 2, function(x) {
  KL.empirical(x, rowSums(barcode.cluster.counts))*sum(x)
})

pts.label <- rep(T, length(cluster.barcode.div))
pdf("dm-ter-bc_barcode_teratoma_bias_correlation.pdf", width = 6, height = 6)
PlotCorrelation(cluster.barcode.div, cluster.div[names(cluster.barcode.div)], box = F, show.corr = F,
                use.label = T, pts.label = pts.label, pt.size = 3, label.font.size = 3,
                pt.color = layer.mapping[names(cluster.barcode.div)]) + 
  theme(legend.position = "none")
dev.off()

pdf("dm-ter-bc_barcode_teratoma_bias_correlation_nolabels.pdf", width = 6, height = 6)
PlotCorrelation(cluster.barcode.div, cluster.div[names(cluster.barcode.div)], box = F, show.corr = F,
                use.label = F, pts.label = pts.label, pt.size = 3.5, label.font.size = 0,
                pt.color = layer.mapping[names(cluster.barcode.div)]) +
  theme(legend.position = "none")
dev.off()


cluster.barcode.div <- sort(cluster.barcode.div)
cluster.barcode.gg <- ggBarplot(cluster.unique.barcodes[names(cluster.barcode.div)], 
                                fill.color = layer.mapping[names(cluster.barcode.div)]) + 
  theme(legend.position = "none") + coord_flip()

cluster.barcode.div.gg <- ggBarplot(cluster.barcode.div, 
                                    fill.color = layer.mapping[names(cluster.barcode.div)]) + 
  theme(axis.text.y = element_blank(), legend.title = element_blank()) +
  coord_flip()

cluster.barcode.div.legend <- get_legend(cluster.barcode.div.gg)
cluster.barcode.div.gg <- cluster.barcode.div.gg + theme(legend.position = "none")

pdf("dm-ter-bc_cluster_unique_barcode_count.pdf", width = 6, height = 6)
plot_grid(cluster.barcode.gg, cluster.barcode.div.gg, align = "h", rel_widths = c(0.7, 0.35))
dev.off()

pdf("dm-ter-bc_cluster_unique_barcode_count_legend.pdf", width = 2, height = 4)
ggdraw(cluster.barcode.div.legend)
dev.off()


## Save results
save.image(output.file)
