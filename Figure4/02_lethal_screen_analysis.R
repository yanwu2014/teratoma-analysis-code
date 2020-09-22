library(swne)
library(perturbLM)
library(Seurat)
library(Matrix)
library(glmnet)
library(ggplot2)

## Output RData file
output.file <- "dm-ter-screen_analysis.RData"


#### Load clustering data ####
seurat.file <- "dm-ter-screen_clustering_Seurat.Robj"
ter <- readRDS(seurat.file)

ter.design.mat <- as.matrix(DesignMatrixGenotypes(UnflattenGroups(ter$teratoma)))
dim(ter.design.mat)

idents <- ter$ident
table(idents)

layer.mapping <- c("MSC/Fib"="Mesoderm",
                   "Retinal Epi"="Ectoderm",
                   "Neurons"="Ectoderm",
                   "Neural Prog"="Ectoderm",
                   "Mid/Hindgut Epi"="Endoderm",
                   "Muscle"="Mesoderm",
                   "Smooth Muscle"="Mesoderm",
                   "Pericytes"="Mesoderm",
                   "Foregut Epi"="Endoderm",
                   "SCP"="Ectoderm",
                   "Kidney Prog"="Mesoderm",
                   "Hematopoietic"="Mesoderm")
layer.mapping <- sort(layer.mapping)
layer.ident <- factor(plyr::revalue(idents, replace = layer.mapping))
table(layer.ident)

## Make UMAP plot
umap.emb <- Embeddings(ter, "umap")

set.seed(135235)
cell.sample <- sample(rownames(umap.emb), size = 20000)
pdf("dm-ter-screen_proj_umap.pdf", width = 5, height = 5)
PlotDims(umap.emb[cell.sample,], sample.groups = idents[cell.sample], 
         show.axes = F, do.label = T, show.legend = F, pt.size = 0.5, 
         alpha = 0.5, seed = 231345)
dev.off()

pdf("dm-ter-screen_proj_umap_nolabels.pdf", width = 5, height = 5)
PlotDims(umap.emb[cell.sample,], sample.groups = idents[cell.sample], 
         show.axes = F, do.label = T, show.legend = F, pt.size = 0.5, 
         alpha = 0.5, seed = 231345, label.size = 0)
dev.off()



#### Run EMD analysis for original screen ####

## Load guides
screen.guide.dict <- ReadGroups("dm-ter-screen-merged_pheno_dict.csv")
screen.guide.dict <- lapply(screen.guide.dict, function(x) intersect(x, colnames(ter)))
names(screen.guide.dict) <- gsub("NKX2-1", "NKX2.1", names(screen.guide.dict))
names(screen.guide.dict) <- gsub("NKX2-5", "NKX2.5", names(screen.guide.dict))
sapply(screen.guide.dict, length)


## Remove any NTC cells that also got a true guide
screen.ko.cells <- unique(unlist(screen.guide.dict[!grepl("NTC", names(screen.guide.dict))], F, F))
screen.ntc.guide.dict <- screen.guide.dict[grepl("NTC", names(screen.guide.dict))]
screen.ntc.guide.dict <- lapply(screen.ntc.guide.dict, function(x) x[!x %in% screen.ko.cells])

# ## Filter for cells that only got one guide
# screen.single.mat <- DesignMatrixGenotypes(screen.guide.dict, max.genotypes = 1, min.cells = 1)
# screen.guide.dict <- lapply(colnames(screen.single.mat), function(g) {
#   g.ix <- screen.single.mat[,g]
#   names(g.ix[g.ix > 0])
# })
# names(screen.guide.dict) <- colnames(screen.single.mat)
# 
# screen.guide.dict <- c(screen.guide.dict, screen.ntc.guide.dict)
# screen.guide.dict[["NTC"]] <- unique(unlist(screen.ntc.guide.dict, F, F))
# sapply(screen.guide.dict, length)
screen.guide.dict[names(screen.ntc.guide.dict)] <- screen.ntc.guide.dict
screen.guide.dict[["NTC"]] <- unique(unlist(screen.ntc.guide.dict, F, F))
sapply(screen.guide.dict, length)

# editing.df <- read.table("dm-ter-screen_gRNA_editing_rates.txt", sep = "\t", header = T)
# editing.df <- subset(editing.df, Editing_Rate > 0.6 & Indel_Rate > 0.5)
# guides.use <- c(as.character(editing.df$Guide), grep("NTC", names(screen.guide.dict), value = T))
guides.use <- c("ASCL1-2",
                "CDX2-1",
                "GATA6-2",
                "HES1-1",
                "GATA3-1",
                "KLF2-1",
                "KLF6-1",
                "MYOG-1",
                "NEUROG1-2",
                "NKX2.1-2",
                "NKX2.5-2",
                "PDX1-2",
                "RUNX1-2",
                "SOX9-1",
                "TULP3-2",
                "TWIST1-2")
guides.use <- c(guides.use, grep("NTC", names(screen.guide.dict), value = T))
screen.guide.dict <- screen.guide.dict[names(screen.guide.dict) %in% guides.use]
sapply(screen.guide.dict, length)


## Make gene dictionary
genes <- unique(sapply(names(screen.guide.dict), ExtractField, field = 1, delim = "-"))
screen.gene.dict <- lapply(genes, function(g) {
  unique(unlist(screen.guide.dict[grepl(g, names(screen.guide.dict))]))
})
names(screen.gene.dict) <- genes

# ## Downsample NTC guides
# set.seed(2347890)
# median.cells <- median(sapply(ter.gene.dict, length))
# ter.gene.dict[ntc.guides] <- lapply(ter.gene.dict[ntc.guides], function(x) sample(x, median.cells))
sapply(screen.gene.dict, length)

## Compute cluster distances
umap.emb <- Embeddings(ter, "umap")
cl.emb <- apply(t(umap.emb), 1, function(x) tapply(x, idents, mean))
cl.dist <- as.matrix(dist(cl.emb))

## Compute gene level EMD
screen.emd <- computeEMD(screen.gene.dict, idents, cl.dist)
screen.emd.ntc <- sort(screen.emd["NTC",])
ggBarplot(screen.emd.ntc, fill.color = "skyblue")


#### Run EMD analysis for screen repool ####

## Load guides
repool.guide.dict <- ReadGroups("dm-ter-screen-repool_pheno_dict_offset.csv")
repool.guide.dict <- lapply(repool.guide.dict, function(x) intersect(x, colnames(ter)))
names(repool.guide.dict) <- gsub("NKX2-1", "NKX2.1", names(repool.guide.dict))
names(repool.guide.dict) <- gsub("NKX2-5", "NKX2.5", names(repool.guide.dict))
sapply(repool.guide.dict, length)

## Remove any NTC cells that also got a true guide
repool.ko.cells <- unique(unlist(repool.guide.dict[!grepl("NTC", names(repool.guide.dict))], F, F))
repool.ntc.guide.dict <- repool.guide.dict[grepl("NTC", names(repool.guide.dict))]
repool.ntc.guide.dict <- lapply(repool.ntc.guide.dict, function(x) x[!x %in% repool.ko.cells])

# ## Filter for cells that only got one guide
# repool.single.mat <- DesignMatrixGenotypes(repool.guide.dict, max.genotypes = 1, min.cells = 1)
# repool.guide.dict <- lapply(colnames(repool.single.mat), function(g) {
#   g.ix <- repool.single.mat[,g]
#   names(g.ix[g.ix > 0])
# })
# names(repool.guide.dict) <- colnames(repool.single.mat)
# 
# repool.guide.dict <- c(repool.guide.dict, repool.ntc.guide.dict)
# repool.guide.dict[["NTC"]] <- unique(unlist(repool.ntc.guide.dict, F, F))
# sapply(repool.guide.dict, length)
repool.guide.dict[names(repool.ntc.guide.dict)] <- repool.ntc.guide.dict
repool.guide.dict[["NTC"]] <- unique(unlist(repool.ntc.guide.dict, F, F))
sapply(repool.guide.dict, length)

## Make gene dictionary
genes <- unique(sapply(names(repool.guide.dict), ExtractField, field = 1, delim = "-"))
repool.gene.dict <- lapply(genes, function(g) {
  unique(unlist(repool.guide.dict[grepl(g, names(repool.guide.dict))]))
})
names(repool.gene.dict) <- genes

# ## Downsample NTC guides
# set.seed(2347890)
# median.cells <- median(sapply(ter.gene.dict, length))
# ter.gene.dict[ntc.guides] <- lapply(ter.gene.dict[ntc.guides], function(x) sample(x, median.cells))
sapply(repool.gene.dict, length)

## Compute cluster distances
umap.emb <- Embeddings(ter, "umap")
cl.emb <- apply(t(umap.emb), 1, function(x) tapply(x, idents, mean))
cl.dist <- as.matrix(dist(cl.emb))

## Compute gene level EMD
repool.emd <- computeEMD(repool.gene.dict, idents, cl.dist)
repool.emd.ntc <- sort(repool.emd["NTC",])
ggBarplot(repool.emd.ntc, fill.color = "skyblue")



#### Run LM analysis with separate screen/repool genes ####

## Create combined screen dictionary
names(screen.gene.dict) <- paste0(names(screen.gene.dict), "_screen")
names(repool.gene.dict) <- paste0(names(repool.gene.dict), "_repool")
combined.gene.dict <- c(screen.gene.dict, repool.gene.dict)
sapply(combined.gene.dict, length)

min.cells <- 20
combined.design.mat <- as.matrix(DesignMatrixGenotypes(combined.gene.dict, max.genotypes = 3, 
                                                       min.cells = 1))
combined.design.mat <- combined.design.mat[,!grepl(":", colnames(combined.design.mat))]
combined.design.mat <- combined.design.mat[,!colnames(combined.design.mat) %in% c("NTC_screen", "NTC_repool")]
combined.design.mat <- combined.design.mat[,colSums(combined.design.mat) >= min.cells]
combined.design.mat <- combined.design.mat[rowSums(combined.design.mat) > 0,]
ter.design.mat <- ter.design.mat[rownames(combined.design.mat),]

combined.cl.cv <- cv.glmnet(cbind(combined.design.mat, ter.design.mat),
                            idents[rownames(combined.design.mat)],
                            alpha = 0, nfolds = 4, nlambda = 40,
                            lambda.min.ratio = 1e-04,
                            family = "multinomial")
plot(combined.cl.cv)

combined.lambda.use <- combined.cl.cv$lambda.1se; combined.lambda.use;
combined.cl.cfs <- CalcClusterEnrich(combined.design.mat, cov.mat = ter.design.mat,
                                     clusters = idents[rownames(combined.design.mat)],
                                     alpha = 0, lambda.use = combined.lambda.use,
                                     ctrl = NULL)


#### Run merged LM analysis ####

genes.cfs.use <- intersect(sapply(names(screen.gene.dict), ExtractField), 
                           sapply(names(repool.gene.dict), ExtractField))
genes.cfs.use <- genes.cfs.use[!genes.cfs.use == "NTC"]

gene.dict <- lapply(genes.cfs.use, function(g) {
  unlist(combined.gene.dict[grepl(g, names(combined.gene.dict))], F, F)
}); names(gene.dict) <- genes.cfs.use
sapply(gene.dict, length)

min.cells <- 20
merged.design.mat <- as.matrix(DesignMatrixGenotypes(gene.dict, max.genotypes = 3, 
                                                     min.cells = 1))
merged.design.mat <- merged.design.mat[,!grepl(":", colnames(merged.design.mat))]
merged.design.mat <- merged.design.mat[,colSums(merged.design.mat) >= min.cells]
merged.design.mat <- merged.design.mat[rowSums(merged.design.mat) > 0,]
ter.design.mat <- as.matrix(DesignMatrixGenotypes(UnflattenGroups(ter$teratoma)))
ter.design.mat <- ter.design.mat[rownames(merged.design.mat),]

merged.cl.cv <- cv.glmnet(cbind(merged.design.mat, ter.design.mat),
                          idents[rownames(merged.design.mat)],
                          alpha = 0, nfolds = 5, nlambda = 50,
                          lambda.min.ratio = 2.5e-04,
                          family = "multinomial")
plot(merged.cl.cv)

merged.lambda.use <- merged.cl.cv$lambda.1se; merged.lambda.use;
merged.cl.cfs <- CalcClusterEnrich(merged.design.mat, cov.mat = ter.design.mat,
                                   clusters = idents[rownames(merged.design.mat)],
                                   alpha = 0, lambda.use = merged.lambda.use,
                                   ctrl = NULL)

merged.cl.p <- CalcClusterEnrichPvals(design.mat = merged.design.mat,
                                      cov.mat = ter.design.mat,
                                      clusters = idents[rownames(merged.design.mat)],
                                      alpha = 0, lambda.use = merged.lambda.use,
                                      n.rand = 4000, n.cores = 8, ctrl = NULL)

merged.cl.fdr <- matrix(p.adjust(as.vector(abs(merged.cl.p)), method = "BH"), 
                        nrow = nrow(merged.cl.p))
rownames(merged.cl.fdr) <- rownames(merged.cl.p)
colnames(merged.cl.fdr) <- colnames(merged.cl.p)

max.sig.fdr <- 0.05
merged.cl.cfs.sig <- merged.cl.cfs
merged.cl.cfs.sig[abs(merged.cl.fdr) > max.sig.fdr] <- 0




#### Compare screen and repool results ####

## Compare EMD across screen and repool
pdf("dm-ter-screen_EMD_correlation.pdf", width = 5, height = 5)
genes.use <- intersect(names(screen.emd.ntc), names(repool.emd.ntc))
PlotCorrelation(screen.emd.ntc[genes.use], repool.emd.ntc[genes.use], box = F, use.label = T,
                pts.label = rep(T, length(genes.use)), show.corr = T, pt.size = 2.25, 
                label.font.size = 5.5, x.lab = "Screen EMD", 
                y.lab = "Repool EMD")
dev.off()


genes.cfs.use <- genes.use[genes.use != "NTC"]
gene.corrs <- sapply(genes.cfs.use, function(g) {
  g1 <- paste0(g, "_screen")
  g2 <- paste0(g, "_repool")
  cor(combined.cl.cfs[,g1], combined.cl.cfs[,g2])
})
sort(gene.corrs, decreasing = T)

## Normalize EMD relative to NTC
gene.avg.emd <- rowMeans(cbind(screen.emd.ntc[guides.use], repool.emd.ntc[guides.use]))
gene.avg.emd <- gene.avg.emd[names(gene.avg.emd) != "NTC"]
emd.scale <- 1/mean(gene.avg.emd[grepl("NTC", names(gene.avg.emd))])
gene.avg.emd <- gene.avg.emd*emd.scale

pdf("dm-ter-screen_emd_vs_reproducibility.pdf", width = 5.5, height = 5)
# saveRDS(gene.avg.emd, file = "dm-ter-screen_avg_emd.Robj")
PlotCorrelation(gene.corrs, gene.avg.emd[names(gene.corrs)], box = F, use.label = T,
                pts.label = rep(T, length(gene.corrs)), show.corr = T, pt.size = 2.25, 
                label.font.size = 5.5, x.lab = "Reproducibility (Correlation)", 
                y.lab = "Avg Effect Size (Earth Mover's Distance)")
dev.off()


sig.genes <- c("TWIST1", "RUNX1", "ASCL1", "CDX2", "KLF6",
               "NTC1", "NTC2", "NTC3", "NTC4", "NTC5")
pdf("dm-ter-screen_gene_cluster_heatmap.pdf", width = 5, height = 3.5)
ggHeat(merged.cl.cfs.sig[names(layer.mapping), sig.genes], x.lab.size = 11, y.lab.size = 11)
dev.off()


gene <- "KLF6"
pdf(paste0("dm-ter-screen_", gene, "_replicate_corr.pdf"), width = 3.25, height = 3.25)
PlotCorrelation(combined.cl.cfs[,paste0(gene, "_screen")], 
                combined.cl.cfs[,paste0(gene, "_repool")], 
                box = F, use.label = T, pts.label = rep(T, nrow(combined.cl.cfs)), 
                show.corr = F, 
                pt.size = 2,
                label.font.size = 4, 
                pt.color = layer.mapping,
                title = NULL,
                x.lab = "Screen Effect (linear model)",
                y.lab = "Repool Effect (linear model)") + 
  theme(legend.position = "none", axis.title = element_blank())
dev.off()

## Get legend
library(cowplot)
ggobj <- PlotCorrelation(guide.cl.cfs[,guide.1], guide.cl.cfs[,guide.2],
                         use.label = T, pts.label = pts.label, pt.size = 1.5, label.font.size = 4.5,
                         show.corr = T, title = gene, box = F, pt.color = layer.mapping)
ggobj_leg <- get_legend(ggobj)

pdf("dm-ter-screen-merged_gRNA_corr_legend.pdf", width = 2, height = 4)
ggdraw(ggobj_leg)
dev.off()

## Plot distribution of cells per gene
screen.gene.cells <- sapply(screen.gene.dict, length)
screen.gene.cells <- screen.gene.cells[!grepl("NTC", names(screen.gene.cells))]
median(screen.gene.cells)

repool.gene.cells <- sapply(repool.gene.dict, length)
repool.gene.cells <- repool.gene.cells[!grepl("NTC", names(repool.gene.cells))]
median(repool.gene.cells)

pdf("dm-ter-screen_orig_cells_per_gene.pdf", width = 4, height = 4)
hist(screen.gene.cells, main = "", xlab = "")
dev.off()

pdf("dm-ter-screen_repool_cells_per_gene.pdf", width = 4, height = 4)
hist(repool.gene.cells, main = "", xlab = "")
dev.off()


## Save results
save.image(output.file)
