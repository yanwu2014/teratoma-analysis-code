library(swne)
library(perturbLM)
library(Seurat)
library(Matrix)
library(DESeq2)
library(ggplot2)


## Set working directory to file directory
setwd("teratoma-analysis-code/Figure4")

## Input/output files
seurat.file <- "../../dm-ter-neural_clustering_Seurat.Robj"
output.file <- "dm-ter-neural-screen_DEG_pseudobulk_byTer.RData"


#### Run Pseudobulk ####

## Load Seurat file
ter <- readRDS(seurat.file)
idents <- ter$ident

## Collapse clusters
tissue.mapping <- c("Radial Glia"="Neural Progenitors",
                    "CycProg" = "Neural Progenitors")
idents <- plyr::revalue(idents, replace = tissue.mapping)

# idents.n <- as.numeric(table(idents))
# names(idents.n) <- names(table(idents))
# idents <- droplevels(idents[idents %in% names(idents.n[idents.n > 200])])
table(idents)

ter <- subset(ter, cells = names(idents))
counts <- GetAssayData(ter, slot = "counts", assay = "RNA")[,colnames(ter)]
counts <- counts[!grepl("mm10", rownames(counts)),]
dim(counts)

idents <- idents[colnames(ter)]
ter.ident <- ter$teratoma

row.metadata <- data.frame(gene_short_name = rownames(counts))
rownames(row.metadata) <- rownames(counts)
rm(ter); invisible(gc());

## Load guides
min.cells <- 100
guide.dict <- ReadGroups("dm-ter-neural-screen_pheno_dict_merged.csv")
guide.dict <- lapply(guide.dict, function(x) {
  x <- gsub("\\.", "_", x)
  intersect(x, colnames(counts))
})
guide.dict <- guide.dict[sapply(guide.dict, length) > min.cells]
sapply(guide.dict, length)

## Remove any NTC cells that also got a true guide
ko.cells <- unique(unlist(guide.dict[!grepl("NTC", names(guide.dict))], F, F))
ntc.guide.dict <- guide.dict[grepl("NTC", names(guide.dict))]
ntc.guide.dict <- lapply(ntc.guide.dict, function(x) x[!x %in% ko.cells])
guide.dict[names(ntc.guide.dict)] <- ntc.guide.dict
sapply(guide.dict, length)


## Global parameters for regression analysis
cl.use <- levels(factor(idents))
print(cl.use)


## Run differential expression analysis
diff.exp.list.byTer <- lapply(cl.use, function(cl) {
  cl.cells <- names(idents[idents == cl])
  cl.guide.dict <- lapply(guide.dict, function(x) intersect(x, cl.cells))
  sapply(cl.guide.dict, length)
  
  ## Format pseudobulk matrix
  cl.counts <- do.call(cbind, lapply(cl.guide.dict, function(cells) {
    batch <- ter.ident[cells]
    do.call(cbind, lapply(levels(batch), function(b) {
      b.cells <- names(batch[batch == b])
      rowSums(as.matrix(counts[,b.cells]))
    }))
  }))
  
  ## Format metadata
  guides <- factor(rep(names(cl.guide.dict), each = nlevels(ter.ident)))
  genes <- sapply(as.character(guides), ExtractField, field = 1, delim = "-")
  ter <- rep(levels(ter.ident), times = length(cl.guide.dict))
  genes <- paste(genes, ter, sep = "_")
  genes[grepl("NTC", genes)] <- "NTC"
  genes <- relevel(factor(genes), ref = "NTC")
  ter <- factor(ter)
  cl.coldata <- data.frame(guides, genes, ter)
  
  ## Run DESeq2
  dds <- DESeqDataSetFromMatrix(countData = cl.counts,
                                colData = cl.coldata,
                                design = ~genes + ter)
  # keep <- rowSums(counts(dds)) >= 10
  # dds <- dds[keep,]
  
  dds <- DESeq(dds)
  # resultsNames(dds) # lists the coefficients
  
  ## Pull out results
  genes.use <- levels(genes)
  genes.use <- genes.use[genes.use != "NTC"]
  genes.res.list <- lapply(genes.use, function(g) {
    res.name <- paste0("genes_", g, "_vs_NTC")
    # res <- lfcShrink(dds, coef = res.name, type = "apeglm")
    # res[order(res$pvalue),]
    results(dds, name = res.name)
  })
  names(genes.res.list) <- genes.use
  
  return(genes.res.list)
})
names(diff.exp.list.byTer) <- cl.use


## Save results
save(diff.exp.list.byTer, file = output.file)



#### Analyze results ####

load("dm-ter-neural-screen_DEG_pseudobulk.RData")
load(output.file)

# diff.exp.list1 <- readRDS("dm-ter-neural-screen_DEG_pseudobulk_ter1.Robj")
# diff.exp.list2 <- readRDS("dm-ter-neural-screen_DEG_pseudobulk_ter2.Robj")

cl <- "Neurons"
gene <- "L1CAM"

df <- diff.exp.list[[cl]][[gene]]
df <- subset(df, padj < 0.1)
sig.genes <- unique(rownames(df))
# fdr.sig.genes <- unique(rownames(subset(df, padj < 0.05)))

# df1 <- diff.exp.list1[["Neurons"]][["L1CAM"]]
# df2 <- diff.exp.list2[["Neurons"]][["L1CAM"]]
df1 <- diff.exp.list.byTer[[cl]][[paste0(gene, "_ter1")]]
df2 <- diff.exp.list.byTer[[cl]][[paste0(gene, "_ter2")]]

lfc1 <- df1$log2FoldChange
names(lfc1) <- rownames(df1)
lfc1 <- lfc1[!is.na(lfc1)]

lfc2 <- df2$log2FoldChange
names(lfc2) <- rownames(df2)
lfc2 <- lfc2[!is.na(lfc2)]

genes.use <- intersect(names(lfc1), names(lfc2))
genes.use <- intersect(genes.use, sig.genes)

lq <- -log(df$padj); names(lq) <- rownames(df);
pts.label <- genes.use %in% rownames(head(df, n = 10))
pdf(paste0("dm-ter-neural_", cl, "_", gene, "_logfc_corr.pdf"), width = 5, height = 4)
PlotCorrelation(lfc1[genes.use], lfc2[genes.use], use.label = T, show.corr = T, pt.color = lq,
                box = F, pts.label = pts.label, pt.size = 1.5, title = paste(gene, cl, sep = "/"),
                x.lab = "Teratoma 1 logFC", y.lab = "Teratoma 2 logFC")
dev.off()
