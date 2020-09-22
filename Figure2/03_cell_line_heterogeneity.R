library(swne)
library(perturbLM)


## Set working directory to file directory
setwd("teratoma-analysis-code/Figure2")


## Assess how well distributed each cell type is across each batch
NormEntropy <- function(cluster.batch.counts) {
  require(entropy)
  
  batch.counts <- colSums(cluster.batch.counts)
  total.counts <- sum(as.numeric(cluster.batch.counts))
  n.batch <- ncol(cluster.batch.counts)
  
  raw.entropy <- sum(apply(cluster.batch.counts, 1, function(x) {
    sum(x)*KL.empirical(x, batch.counts)
  }))
  
  return(1 - raw.entropy/(log(n.batch)*total.counts))
}


## Load data
h1.meta <- read.table("../Figure1/dm-ter-human-merged_metadata.tsv", header = T, sep = "\t")
cell.line.meta <- read.table("dm-ter-chimera-merged_metadata.tsv", header = T, sep = "\t")
organoid.meta <- read.table("../Reference_Data/velasco_organoid_meta_combined.txt", header = T, sep = "\t")


## Quantify H1 cell heterogeneity using MI
h1.clusters <- h1.meta$cluster; names(h1.clusters) <- h1.meta$cell;
h1.clusters.list <- UnflattenGroups(h1.clusters)

h1.batch <- h1.meta$teratoma; names(h1.batch) <- h1.meta$cell;
h1.batch.list <- UnflattenGroups(h1.batch)

h1.ter.cluster.counts <- GroupOverlapCounts(h1.clusters.list, h1.batch.list)
h1.norm.entropy <- NormEntropy(h1.ter.cluster.counts)


## Quantify cell line heterogeneity using MI
cell.line.clusters <- cell.line.meta$ident; names(cell.line.clusters) <- cell.line.meta$cell;
cell.line.clusters.list <- UnflattenGroups(cell.line.clusters)

cell.line.batch <- cell.line.meta$cell_line; names(cell.line.batch) <- cell.line.meta$cell;
cell.line.batch.list <- UnflattenGroups(cell.line.batch)

line.ter.cluster.counts <- GroupOverlapCounts(cell.line.clusters.list, cell.line.batch.list)
line.norm.entropy <- NormEntropy(line.ter.cluster.counts)


## Quantify organoid heterogeneity using MI
organoid.meta <- subset(organoid.meta, Batch == "PGP1_6mon_Batch1")

organoid.clusters <- organoid.meta$CellType; names(organoid.clusters) <- organoid.meta$NAME;
organoid.clusters.list <- UnflattenGroups(organoid.clusters)

organoid.batch <- factor(organoid.meta$Organoid); names(organoid.batch) <- organoid.meta$NAME;
levels(organoid.batch) <- c("org16", "org17", "org18")
organoid.batch.list <- UnflattenGroups(organoid.batch)

organoid.ter.cluster.counts <- GroupOverlapCounts(organoid.clusters.list, organoid.batch.list)
organoid.norm.entropy <- NormEntropy(organoid.ter.cluster.counts)

norm.entropy <- c(h1.norm.entropy, line.norm.entropy, organoid.norm.entropy)
names(norm.entropy) <- c("H1 Teratomas", "Cell Line Teratomas", "Brain Organoids")

pdf("teratoma_organoid_heterogeneity_barplot.pdf", width = 2, height = 4)
ggBarplot(norm.entropy)
dev.off()


