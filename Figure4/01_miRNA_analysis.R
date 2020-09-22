library(Seurat)
library(swne)
library(perturbLM)
library(ggplot2)


## Set working directory to file directory
setwd("teratoma-analysis-code/Figure4")

## miR-124 Analysis ##

## Output file
output.file <- "dm-ter-miRNA_analysis.RData"


## Load datasets
counts.dirs <- c("../Counts/miR-124-matrices/dm-ter-miRNA-pos-it1", 
                 "../Counts/miR-124-matrices/dm-ter-miRNA-pos-it2", 
                 "../Counts/miR-124-matrices/dm-ter-miRNA-pos-ipit1", 
                 "../Counts/miR-124-matrices/dm-ter-miRNA-pos-ipit2",
                 "../Counts/miR-124-matrices/dm-ter-miRNA-neg")
class.files <- c("../Counts/miR-124-matrices/dm-ter-miRNA-pos-it1_species_class.csv", 
                 "../Counts/miR-124-matrices/dm-ter-miRNA-pos-it2_species_class.csv", 
                 "../Counts/miR-124-matrices/dm-ter-miRNA-pos-ipit1_species_class.csv", 
                 "../Counts/miR-124-matrices/dm-ter-miRNA-pos-ipit2_species_class.csv",
                 "../Counts/miR-124-matrices/dm-ter-miRNA-neg_species_class.csv")

counts.list <- lapply(1:length(counts.dirs), function(i) {
  Read10XHuman(counts.dirs[[i]], class.files[[i]])
})
names(counts.list) <- c("IT1", "IT2", "IPIT1", "IPIT2", "Neg")

## Map teratoma screen cells to reference clusters
ref.seurat.file <- "../Figure1/dm-ter-human-merged_clustering_v3.seurat.Robj"

## Load data
ref <- readRDS(ref.seurat.file)
DefaultAssay(ref) <- "integrated"

## Load metadata
meta.data <- read.table("../Figure1/dm-ter-human-merged_metadata.tsv", header = T, sep = "\t")
ref.ident <- meta.data$cluster; names(ref.ident) <- meta.data$cell;
table(ref.ident)

## Collapse clusters
tissue.mapping <- c("Airway Epi"="Foregut Epi",
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
tissue.ident <- plyr::revalue(ref.ident, replace = tissue.mapping)
table(tissue.ident)

## Create Seurat objects and run label transfer
obj.list <- lapply(counts.list, function(counts) {
  ter <- CreateSeuratObject(counts, min.cells = 10, min.features = 200)
  ter <- NormalizeData(ter)
  ter <- ScaleData(ter)
  
  transfer.anchors <- FindTransferAnchors(reference = ref, ter, dims = 1:20)
  chimera.ident.df <- TransferData(anchorset = transfer.anchors, refdata = tissue.ident, dims = 1:20)
  ter$ident <- chimera.ident.df$predicted.id
  
  return(ter)
})

## Order idents by germ layer
tissue.ident.order <- c("Early Neurons", "Neuronal Prog", "SCP", "Retinal Epi",
                        "MSC/Fib", "Muscle", "Hematopoietic", "Smooth Muscle", "Pericytes",
                        "Kidney Prog", "Foregut Epi", "Mid/Hindgut Epi")

## Make table of cluster by teratoma cell counts
ter.ident.list <- lapply(obj.list, function(ter) {
  x <- rep(0, length(table(tissue.ident)))
  names(x) <- names(table(tissue.ident))
  x[names(table(ter$ident))] <- table(ter$ident)
  x
})
ter.ident.tbl <- do.call(cbind, ter.ident.list)
ter.ident.frac <- t(t(ter.ident.tbl)/colSums(ter.ident.tbl))

ter.logfc <- sapply(c("IT1", "IT2", "IPIT1", "IPIT2"), function(id) {
  log2((ter.ident.frac[,id] + 1e-3)/(ter.ident.frac[,"Neg"] + 1e-3))
})

# pdf("dm-ter-miRNA_logfc_heatmap.pdf", width = 3.75, height = 4.25)
heatmap.gg <- ggHeat(ter.logfc[tissue.ident.order,], x.lab.size = 11, y.lab.size = 11) + 
  theme(legend.position = "none")
# dev.off()

library(ggpubr)
heatmap.leg <- ggpubr::get_legend(ggHeat(ter.logfc[tissue.ident.order,], x.lab.size = 11, y.lab.size = 11))

## Build null distribution with H1 teratoma cells to get variance
tissue.ident.list <- UnflattenGroups(tissue.ident)
ter.ident.list <- UnflattenGroups(ref$batch)

ref.ident.tbl <- GenotypeClusterCounts(tissue.ident.list, ter.ident.list)
ref.ident.frac <- t(t(ref.ident.tbl)/colSums(ref.ident.tbl))
ref.ident.var <- apply(ref.ident.frac, 1, function(x) var(x))

## Compute z-scores
ter.it.avg.frac <- rowMeans(ter.ident.frac[,c("IT1", "IT2")])
ter.it.var <- apply(ter.ident.frac[,c("IT1", "IT2")], 1, var)

## Number of teratomas used to compute pooled standard deviation
num.ref.ter <- ncol(ref.ident.tbl)
num.pos.ter <- 2 ## 2 IT and IPIT teratomas each

## Use Cohen's pooled standard deviation
ter.it.pooled.sd <- sqrt((ter.it.var + (num.ref.ter - 1)*ref.ident.var)/num.ref.ter)

ter.it.z <- (ter.it.avg.frac - ter.ident.frac[,"Neg"])/ter.it.pooled.sd
ter.it.z[ter.it.z < -6] <- -6

# pdf("dm-ter-miRNA_zscore_barplot_IT.pdf", width = 4, height = 4.25)
ggBarplot(ter.gcv.z[tissue.ident.order], fill.color = "skyblue") + coord_flip()
# dev.off()

ter.ipit.avg.frac <- rowMeans(ter.ident.frac[,c("IPIT1", "IPIT2")])
ter.ipit.var <- apply(ter.ident.frac[,c("IPIT1", "IPIT2")], 1, var)

## Cohen's pooled standard deviation
ter.ipit.pooled.sd <- sqrt((ter.ipit.var + (num.ref.ter - 1)*ref.ident.var)/num.ref.ter)

ter.ipit.z <- (ter.ipit.avg.frac - ter.ident.frac[,"Neg"])/ter.ipit.pooled.sd
ter.ipit.z[ter.ipit.z < -6] <- -6

# pdf("dm-ter-miRNA_zscore_barplot_IPIT.pdf", width = 4.5, height = 4.5)
ggBarplot(ter.ipit.gcv.z[tissue.ident.order], fill.color = "skyblue") + coord_flip()
# dev.off()

ter.z <- rbind(data.frame(z = ter.it.z, cluster = names(ter.it.z), inject = "IT"),
               data.frame(z = ter.ipit.z, cluster = names(ter.ipit.z), inject = "IPIT"))
ter.z$cluster <- factor(ter.z$cluster, levels = tissue.ident.order)


barplot.gg <- ggplot(ter.z, aes(x = cluster, y = z, fill = inject)) +
  geom_bar(stat = 'identity', position = 'dodge') + 
  theme_classic() + 
  theme(axis.title = element_blank(), legend.title = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_text(size = 12),
        legend.position = "none") + 
  coord_flip()

barplot.leg <- ggpubr::get_legend(ggplot(ter.z, aes(x = cluster, y = z, fill = inject)) +
                                    geom_bar(stat = 'identity', position = 'dodge') + 
                                    theme_classic() + 
                                    theme(legend.text = element_text(size = 11)) + 
                                    coord_flip())

# pdf("dm-ter-miRNA_heatmap_barplot.pdf", width = 5, height = 4)
library(cowplot)
plot_grid(heatmap.gg, barplot.gg, align = "h", rel_widths = c(0.6,0.4))
# dev.off()


# pdf("dm-ter-miRNA_heatmap_legend.pdf", width = 2, height = 4)
as_ggplot(heatmap.leg)
# dev.off()

# pdf("dm-ter-miRNA_barplot_legend.pdf", width = 2, height = 4)
as_ggplot(barplot.leg)
# dev.off()


## Save output
save.image(output.file)

