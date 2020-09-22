library(Seurat)
library(swne)
library(cellMapper)

#### Prep MCA dataset ####

load("Reference_Data/mca_full_debatched.RData")
metadata <- read.table("Reference_Data/mca_cell_annotation.csv", sep = ",", header = T)
rownames(metadata) <- metadata$Cell.name

mca.counts <- mca.counts[,colnames(mca.counts) %in% rownames(metadata)]
metadata <- droplevels(metadata[colnames(mca.counts),])
colnames(metadata) <- c("cell_name", "clusterID", "tissue", "batch", "cell_barcode", "cell_type", "mouse_id")

## Map mouse genes to human genes
mouse2human <- MouseHumanMapping(rownames(mca.counts))
mca.counts <- mca.counts[names(mouse2human),]; rownames(mca.counts) <- mouse2human[rownames(mca.counts)];
mca.counts <- mca.counts[!duplicated(rownames(mca.counts)),]

## Filter reference counts
mca.counts <- FilterData(mca.counts, min.samples.frac = 1e-4, trim = 1e-5, min.nonzero.features = 200)
dim(mca.counts)

## Create Seurat object
mca <- CreateSeuratObject(mca.counts, project = "MouseCellAtlas", normalization.method = "LogNormalize",
                              scale.factor = median(Matrix::colSums(mca.counts)))
mca <- AddMetaData(mca, metadata = metadata[colnames(mca.counts),])
mca <- SetAllIdent(mca, id = "cell_type")

## Save reference classifier
# save(mca, file = "Reference_Data/mca_seurat_mouse.RData")
save(mca, file = "Reference_Data/mca_seurat_human.RData")


#### Prep MOCA dataset ####

## Load data
cds <- readRDS("Reference_Data/cds_cleaned_sampled_100k.RDS")
markers.df <- read.table("Reference_Data/DE_gene_main_cluster.csv", header = T, sep = ",")

## Extract counts, metadata, and gene metadata
counts <- cds@assayData$exprs
metadata <- cds@phenoData@data
gene.data <- cds@featureData@data 
rownames(counts) <- as.character(gene.data$gene_short_name)

## Map mouse genes to human genes
mouse2human <- MouseHumanMapping(rownames(counts))
counts <- counts[names(mouse2human),]
rownames(counts) <- mouse2human[rownames(counts)]

markers.df <- subset(markers.df, qval < 1e-3 & fold.change > 1.5)
all.mouse.markers <- as.character(unique(markers.df$gene_short_name))
all.mouse.markers <- intersect(all.mouse.markers, names(mouse2human))
all.human.markers <- mouse2human[all.mouse.markers]

## Create Seurat object
moca <- CreateSeuratObject(counts, project = "MOCA", meta.data = metadata)

## Save reference classifier
save(moca, file = "Reference_Data/moca_seurat_human.RData")
