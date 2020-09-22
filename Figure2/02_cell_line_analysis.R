library(methods)
library(Seurat)
library(swne)
library(ggplot2)
library(perturbLM)
library(cowplot)
library(entropy)


## Set working directory to file directory
setwd("teratoma-analysis-code/Figure2")

## Output file
output.file <- "dm-ter-human_cell_line_analysis.RData"


#### Make teratoma cell type stacked barplot ####
chimera.metadata <- read.table("dm-ter-chimera-merged_metadata.tsv", header = T, sep = "\t")
cell.line <- chimera.metadata$cell_line; names(cell.line) <- chimera.metadata$cell;
chimera.ident <- chimera.metadata$ident; names(chimera.ident) <- chimera.metadata$cell;

line.cluster.counts <- GroupOverlapCounts(UnflattenGroups(chimera.ident), 
                                          UnflattenGroups(cell.line))


## Load cell type and teratoma assignments
meta.data <- read.table("../Figure1/dm-ter-human-merged_metadata.tsv", header = T, 
                        row.names = NULL, sep = "\t")

h1.clusters <- meta.data$cluster; names(h1.clusters) <- meta.data$cell;
h1.clusters.list <- UnflattenGroups(h1.clusters)

h1.ter.id <- meta.data$teratoma; names(h1.ter.id) <- meta.data$cell;
h1.ter.id.list <- UnflattenGroups(h1.ter.id)

h1.cluster.counts <- GroupOverlapCounts(h1.clusters.list, h1.ter.id.list)
colnames(h1.cluster.counts) <- paste0("H1-", colnames(h1.cluster.counts))

ter.batch.counts <- t(cbind(h1.cluster.counts, line.cluster.counts[rownames(h1.cluster.counts),]))
ter.batch.counts.frac <- ter.batch.counts/rowSums(ter.batch.counts)
ter.cluster.df <- setNames(reshape2::melt(ter.batch.counts.frac), c("Teratoma", "Cluster", "Frac"))

pdf("dm-ter-human-merged_teratoma_celltype_barplot.pdf", width = 9, height = 5)
ggplot(ter.cluster.df) +
  geom_bar(aes(x = Teratoma, y = Frac, fill = Cluster), stat = "identity", 
           width = 0.75, color = "black") + 
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 14, angle = 90, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 13), legend.title = element_blank())
dev.off()


msc.frac <- rowSums(ter.batch.counts.frac[,grepl("MSC", colnames(ter.batch.counts))])

pdf("dm-ter-human-merged_MSC_fraction_barplot.pdf", width = 6, height = 3)
ggBarplot(msc.frac, y.lim = c(0,1))
dev.off()

write.table(t(ter.batch.counts), file = "dm-ter-human-merged_ter_cluster_counts.tsv",
            sep = "\t", quote = F)
write.table(t(ter.batch.counts.frac), file = "dm-ter-human-merged_ter_cluster_frac.tsv",
            sep = "\t", quote = F)



#### Make Teratoma Germ Layer Barplot ####

## Load cell type to germ layer mappings
layer.mapping.df <- read.table("../Figure1/dm-ter-human-merged_cluster_layer_mapping.txt", header = T, sep = "\t")
layer.mapping <- as.character(layer.mapping.df$Germ_Layer); names(layer.mapping) <- layer.mapping.df$Cell_Type;


zebrafish <- c(59817, 1375, 20116) ## Output from the 01_count_zebrafish_cells.R script in Reference_Data
names(zebrafish) <- c("Ectoderm", "Endoderm", "Mesoderm")


moca.anno <- read.table("../Reference_Data/moca_mapped_cell_annotations.tsv", header = T, sep = "\t")
mouse <- as.numeric(table(moca.anno$Germ_Layer))
names(mouse) <- names(table(moca.anno$Germ_Layer))
mouse <- mouse[c("Ectoderm", "Endoderm", "Mesoderm")]


library(data.table)
chimera.layer <- plyr::revalue(chimera.ident, replace = layer.mapping)
line.layer.counts <- GroupOverlapCounts(UnflattenGroups(chimera.layer), 
                                        UnflattenGroups(cell.line))

h1.layer <- plyr::revalue(h1.clusters, replace = layer.mapping)
h1.layer.counts <- GroupOverlapCounts(UnflattenGroups(h1.layer), h1.ter.id.list)
colnames(h1.layer.counts) <- paste0("H1-", colnames(h1.layer.counts))

ter.layer.counts <- t(cbind(h1.layer.counts, line.layer.counts[rownames(h1.layer.counts),]))
ter.layer.counts <- rbind(ter.layer.counts, zebrafish[colnames(ter.layer.counts)], 
                          mouse[colnames(ter.layer.counts)])
rownames(ter.layer.counts)[[11]] <- "Zebrafish Embryo"
rownames(ter.layer.counts)[[12]] <- "Mouse Embryo"

ter.layer.counts.frac <- ter.layer.counts/rowSums(ter.layer.counts)
ter.layer.df <- setNames(reshape2::melt(ter.layer.counts.frac), c("Teratoma", "Layer", "Frac"))

pdf("dm-ter-human-merged_teratoma_layer_barplot.pdf", width = 7, height = 4.5)
ggplot(ter.layer.df) +
  geom_bar(aes(x = Teratoma, y = Frac, fill = Layer), stat = "identity", 
           width = 0.75, color = "black") + 
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 14, angle = 90, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 13), legend.title = element_blank())
dev.off()


#### Make teratoma cell type stacked barplot in Figure S2B ####
ter.cluster.counts <- GroupOverlapCounts(h1.clusters.list, h1.ter.id.list)
ter.cluster.frac <- t(ter.cluster.counts)/colSums(ter.cluster.counts)
ter.cluster.frac <- t(t(ter.cluster.frac)/colSums(ter.cluster.frac))
ter.counts <- colSums(ter.cluster.counts)

cluster.div <- apply(ter.cluster.counts, 1, function(x) {
  KL.empirical(x, colSums(ter.cluster.counts))*sum(x)
})
cluster.div <- sort(cluster.div)
cluster.div.bar <- ggBarplot(cluster.div, fill.color = "lightgrey") + coord_flip() +
  theme(axis.text.y = element_blank())

ter.cluster.frac <- ter.cluster.frac[,names(cluster.div)]
ter.cluster.df <- setNames(reshape2::melt(ter.cluster.frac), c("Teratoma", "Cluster", "Frac"))

stacked.bar <- ggplot(ter.cluster.df) +
  geom_bar(aes(x = Cluster, y = Frac, fill = Teratoma), stat = "identity", 
           width = 0.75, color = "black") + 
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(hjust = 1, size = 14, angle = 90, color = "black"),
        axis.text.y = element_text(size = 16, color = "black")) + 
  scale_fill_brewer(palette = "Set2") + coord_flip()



pdf("dm-ter-human-merged_cluster_teratoma_barplot.pdf", height = 7, width = 8)
plot_grid(stacked.bar + theme(legend.position = "none"), cluster.div.bar, align = "h", 
          rel_widths = c(0.65, 0.35))
dev.off()

## Plot legend separately
legend <- get_legend(stacked.bar)

pdf("dm-ter-human-merged_cluster_teratoma_barplot_legend.pdf", height = 7.5, width = 2.5)
ggdraw(legend)
dev.off()

## Save results
save.image(output.file)

