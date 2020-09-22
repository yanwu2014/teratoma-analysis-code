library(swne)
library(perturbLM)
library(Seurat)
library(Matrix)
library(DESeq2)
library(snow)
library(liger)


## Set working directory to file directory
setwd("teratoma-analysis-code/Figure4")

## Input/output files
output.file <- "dm-ter-neural-screen_DEG_pseudobulk.RData"
seurat.file <- "dm-ter-neural_clustering_Seurat.Robj"


#### Run Pseudobulk ####

## Load Seurat file
ter <- readRDS(seurat.file)
idents <- ter$ident
# idents <- Idents(ter)

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

col.metadata <- data.frame(ident = idents, ter = ter$teratoma, run = ter$run)
rownames(col.metadata) <- colnames(ter)

ter.ident <- col.metadata$ter
names(ter.ident) <- rownames(col.metadata)

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

# guide.dict[names(ntc.guide.dict)] <- NULL
# sapply(guide.dict, length)

## Global parameters for regression analysis
cl.use <- levels(factor(idents))
print(cl.use)

diff.exp.list <- lapply(cl.use, function(cl) {
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
  genes[grepl("NTC", genes)] <- "NTC"
  genes <- relevel(factor(genes), ref = "NTC")
  ter <- factor(rep(levels(ter.ident), times = length(cl.guide.dict)))
  cl.coldata <- data.frame(guides, genes, ter)
  
  ## Run DESeq2
  dds <- DESeqDataSetFromMatrix(countData = cl.counts,
                                colData = cl.coldata,
                                design= ~ ter + genes)
  dds <- DESeq(dds)
  resultsNames(dds) # lists the coefficients
  
  ## Pull out results
  genes.use <- levels(genes)
  genes.use <- genes.use[genes.use != "NTC"]
  genes.res.list <- lapply(genes.use, function(g) {
    res.name <- paste0("genes_", g, "_vs_NTC")
    res <- lfcShrink(dds, coef = res.name, type = "apeglm")
    res[order(res$pvalue),]
  })
  names(genes.res.list) <- genes.use
  
  return(genes.res.list)
})
names(diff.exp.list) <- cl.use


## Save results
save(diff.exp.list, file = output.file)



#### Analyze results ####
load(output.file)

sig.diff.exp.list <- lapply(diff.exp.list, function(guide.list) {
  do.call(rbind,lapply(names(guide.list), function(g) {
    df <- guide.list[[g]]
    df$gene <- rownames(df)
    df$guide <- g
    df <- subset(df, padj < 0.1)
    df[order(df$pvalue),]
  }))
})

genes.use <- unique(unlist(lapply(sig.diff.exp.list, function(df) df$guide), F, F))
genes.use <- genes.use[!grepl("NTC", genes.use)]
sig.diff.exp.n <- do.call(cbind, lapply(sig.diff.exp.list, function(df) {
  n <- rep(0, length(genes.use)); names(n) <- genes.use;
  n[names(table(df$guide))] <- table(df$guide)
  n
}))
sig.diff.exp.n
write.table(sig.diff.exp.n, file = "dm-ter-neural_diff_exp_pseudobulk_summary.tsv", sep = "\t", quote = F)

out.dir <- "DEG_Results/"
dir.create(out.dir)
for(cl in names(sig.diff.exp.list)) {
  cl.out.file <- paste0(out.dir, "dm-ter-neural_diff_exp_pseudobulk_", gsub("/| ", "\\.", cl), "_table.tsv")
  write.table(sig.diff.exp.list[[cl]], file = cl.out.file, sep = "\t", quote = F, row.names = F)
}


## Inspect results
df <- diff.exp.list[["Neurons"]][["L1CAM"]]
df <- df[order(df$pvalue),]
head(df, n = 10)