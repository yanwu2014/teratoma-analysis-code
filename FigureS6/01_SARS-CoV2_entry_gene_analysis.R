library(Seurat)

## Read input
ter <- readRDS("../Figure1/dm-ter-human-merged_clustering_v3.seurat.Robj")
ter.metadata <- read.table("../Figure1/dm-ter-human-merged_metadata.tsv", header = T, sep = "\t")
named.ident <- ter.metadata$cluster; names(named.ident) <- ter.metadata$cell
named.ident <- named.ident[colnames(ter)]


## Visualize expression of COVID entry genes ACE2 and TMPRSS2
covid.genes <- c("ACE2", "TMPRSS2")
covid.gene.counts <- GetAssayData(ter, assay = "RNA")[covid.genes,]
covid.gene.cluster.expr <- apply(covid.gene.counts, 1, function(x) tapply(x, named.ident, mean))
covid.gene.cluster.df <- reshape2::melt(covid.gene.cluster.expr)

pdf("dm-ter-human-merged_covid_gene_expr_barplot.pdf", width = 6, height = 4)
ggplot(data = covid.gene.cluster.df, aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7, color = "black") + 
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 12, angle = 90, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.title = element_blank(), legend.position = "top")
dev.off()


gut.cells <- names(named.ident[named.ident == "Mid/Hindgut Epi"])
ace2.expr <- covid.gene.counts["ACE2",gut.cells]
tmprss2.expr <- covid.gene.counts["TMPRSS2",gut.cells]

library(perturbLM)
pdf("dm-ter-human-merged_covid_gene_expr_mid_hindgut_correlation.pdf", width = 4.5, height = 4)
PlotCorrelation(ace2.expr, tmprss2.expr, use.label = F, x.lab = "Mid/Hindgut ACE2 Expr", 
                y.lab = "Mid/Hindgut TMPRSS2 Expr")
dev.off()