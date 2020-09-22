library(Seurat)
library(swne)
library(cellMapper)
library(ggplot2)
library(perturbLM)
library(ggplot2)

## Set working directory to file directory
setwd("teratoma-analysis-code/Figure1")

## Input/output files
output.file <- "dm-ter-human-merged_mouse_fetal_mapping.RData"
input.file <- "dm-ter-human-merged_clustering_v3.seurat.Robj"

## These reference files were generated from the Mouse Cell Atlas data
## http://bis.zju.edu.cn/MCA/
## These paths may change depending on where you download your data to
ref.file <- "../Reference_Data/moca_seurat_human.RData"
ref.markers.file <- "../Reference_Data/moca_DEGs_main_cluster.csv"

## Load reference data
load(ref.file)
moca <- NormalizeData(moca)

ref.markers.df <- read.table(ref.markers.file, header = T, sep = ",", stringsAsFactors = F)
mouse2human <- MouseHumanMapping(as.character(ref.markers.df$gene_short_name))
ref.markers.df <- subset(ref.markers.df, gene_short_name %in% names(mouse2human))
ref.markers.df$human_gene_name <- mouse2human[ref.markers.df$gene_short_name]
ref.markers.df <- subset(ref.markers.df, qval < 1e-8)
ref.genes.use <- unique(ref.markers.df$human_gene_name)

ref.metadata <- read.table("../Reference_Data/moca_cell_annotate.csv", header = T,
                           sep = ",", stringsAsFactors = F)
ref.metadata <- subset(ref.metadata, sample %in% colnames(moca))
ref.clusters <- factor(ref.metadata$Main_cell_type)
names(ref.clusters) <- ref.metadata$sample

# ref.sub.metadata <- subset(ref.metadata, Main_cell_type == "Endothelial cells")
# table(ref.sub.metadata$Sub_trajectory_name)
# head(ref.sub.metadata)

ref.sub.traj <- ref.metadata$Sub_trajectory_name
names(ref.sub.traj) <- ref.metadata$sample

moca$Main_Cell_Type <- ref.clusters[colnames(moca)]
moca$Sub_Trajectory <- ref.sub.traj[colnames(moca)]
moca <- SetIdent(moca, value = "Main_Cell_Type")
table(moca$Main_Cell_Type)

## Load teratoma data
ter <- readRDS(input.file)
DefaultAssay(ter) <- "RNA"

ter.metadata <- read.table("dm-ter-human-merged_metadata.tsv", header = T, sep = "\t")
rownames(ter.metadata) <- ter.metadata$cell
ter$named.ident <- ter.metadata[colnames(ter), "cluster"]
ter <- SetIdent(ter, value = ter$named.ident)

## Select genes to train on
ter.markers.df <- read.table("dm-ter-human-merged_named_cluster_markers.tsv", header = T, sep = "\t",
                             stringsAsFactors = F)
ter.markers.df <- subset(ter.markers.df, avg_logFC > 0.25)
ter.genes.use <- as.character(ter.markers.df$gene)

# moca <- FindVariableFeatures(moca, nfeatures = 4e3)
# ref.genes.use <- VariableFeatures(moca)
genes.use <- union(ref.genes.use, ter.genes.use)
genes.use <- genes.use[genes.use %in% rownames(moca) & genes.use %in% rownames(ter)]

## Load key markers
key.markers.df <- read.table("dm-ter-human-merged_marker_cluster_mapping.txt", header = T, sep = "\t")
min.markers <- unique(key.markers.df$Marker)
cl.order <- unique(key.markers.df$Cell.Type)

## Scale reference
VariableFeatures(moca) <- genes.use
moca <- ScaleData(moca, features = genes.use)

## Scale query dataset
ter <- ScaleData(ter, features = genes.use)

## Run cluster correlations
correlation.res <- MapClustersCor(GetAssayData(ter, slot = "scale.data"), ter$named.ident,
                                  GetAssayData(moca, slot = "scale.data"), moca$Main_Cell_Type,
                                  genes.use = genes.use, metric = "pearson")

ter.ref.correlations <- correlation.res$cluster.correlations
ter.ref.correlations <- ter.ref.correlations[,cl.order]


pdf("dm-ter-human-merged_moca_tissue_mapping.pdf", width = 6.75, height = 9)
ggHeat(ter.ref.correlations, clustering = "row", x.lab.size = 10, y.lab.size = 10) + 
  theme(legend.position = "top")
dev.off()

## Save output
save.image(output.file)


ter.cluster <- "Foregut"
ref.cluster <- "Epithelial cells"

ter.cl.expr <- correlation.res$query.cluster.data[,ter.cluster]
ref.cl.expr <- correlation.res$ref.cluster.data[,ref.cluster]
pts.label <- ter.cl.expr > 1 & ref.cl.expr > 1
PlotCorrelation(ter.cl.expr, ref.cl.expr, pts.label = pts.label,
                show.corr = T, use.label = T)

