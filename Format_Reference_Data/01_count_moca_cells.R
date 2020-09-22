library(perturbLM)

## Load data
metadata <- read.table("Reference_Data/moca_cell_annotate.csv", sep = ",", header = T)
levels(metadata$Main_trajectory)

main.cluster <- metadata$Main_cell_type
names(main.cluster) <- metadata$sample
levels(main.cluster)


#### Count cells in each germ layer for MOCA ####
layer.mapping <- c("Cardiac muscle lineages" = "Mesoderm",
                   "Cholinergic neurons" = "Ectoderm",
                   "Chondroctye progenitors" = "Mesoderm",
                   "Chondrocytes & osteoblasts" = "Mesoderm",
                   "Connective tissue progenitors" = "Mesoderm",
                   "Definitive erythroid lineage" = "Mesoderm",
                   "Early mesenchyme" = "Mesoderm",
                   "Endothelial cells" = "Mesoderm",
                   "Ependymal cell" = "Ectoderm",
                   "Epithelial cells" = NA,
                   "Excitatory neurons" = "Ectoderm",
                   "Granule neurons" = "Ectoderm",
                   "Hepatocytes" = "Endoderm",
                   "Inhibitory interneurons" = "Ectoderm",
                   "Inhibitory neuron progenitors" = "Ectoderm",
                   "Inhibitory neurons" = "Ectoderm",
                   "Intermediate Mesoderm" = "Mesoderm",
                   "Isthmic organizer cells" = "Ectoderm",
                   "Jaw and tooth progenitors" = "Ectoderm",
                   "Lens" = "Ectoderm",
                   "Limb mesenchyme" = "Mesoderm",
                   "Megakaryocytes" = "Mesoderm",
                   "Melanocytes" = "Ectoderm",
                   "Myocytes" = "Mesoderm",
                   "Neural progenitor cells" = "Ectoderm",
                   "Neural Tube" = "Ectoderm",
                   "Neutrophils" = "Mesoderm",
                   "Notochord cells" = "Ectoderm",
                   "Oligodendrocyte Progenitors" = "Ectoderm",
                   "Osteoblasts" = "Mesoderm",
                   "Postmitotic premature neurons" = "Ectoderm",
                   "Premature oligodendrocyte" = "Ectoderm",
                   "Primitive erythroid lineage" = "Mesoderm",
                   "Radial glia" = "Ectoderm",
                   "Schwann cell precursor" = "Ectoderm",
                   "Sensory neurons" = "Ectoderm",
                   "Stromal cells" = "Mesoderm",
                   "White blood cells" = "Mesoderm")

layer <- as.character(plyr::revalue(main.cluster, replace = layer.mapping))
names(layer) <- metadata$sample
table(layer)

epi.metadata <- subset(metadata, Main_cell_type == "Epithelial cells")
epi.subcluster <- droplevels(epi.metadata$Sub_trajectory_name)
names(epi.subcluster) <- epi.metadata$sample
levels(epi.subcluster)

epi.layer.mapping <- c("Apical ectodermal ridge trajectory" = "Ectoderm",
                       "Auditory epithelial trajectory" = NA,
                       "Branchial arch epithelial trajectory" = NA,
                       "Epidermis trajectory" = "Ectoderm",
                       "Lung epithelial trajectory" = "Endoderm",
                       "Midgut/Hindgut epithelial trajectory" = "Endoderm",
                       "Olfactory epithelial trajectory" = "",
                       "Pericardium trajectory" = "Endoderm",
                       "Primordial germ cell trajectory" = NA,
                       "Renal epithelial trajectory" = "Mesoderm",
                       "Retina epithelial trajectory" = "Ectoderm",
                       "Stomach epithelial trajectory" = "Endoderm",
                       "Urothelium trajectory" = NA)

epi.layer <- plyr::revalue(epi.subcluster, replace = epi.layer.mapping)
names(epi.layer) <- epi.metadata$sample

layer[names(epi.layer)] <- as.character(epi.layer)
table(layer)

# mapped.metadata <- data.frame(Cell = metadata$sample, Cell_Type = main.cluster, Germ_Layer = layer)
# write.table(mapped.metadata, file = "moca_mapped_cell_annotations.tsv", sep = "\t")



#### Count cells in each cluster for MOCA ####
main.cluster.epi <- as.character(main.cluster)
names(main.cluster.epi) <- names(main.cluster)
main.cluster.epi[names(epi.subcluster)] <- as.character(epi.subcluster)
main.cluster.epi <- factor(main.cluster.epi)

main.cluster.epi.counts <- as.numeric(table(main.cluster.epi))
names(main.cluster.epi.counts) <- names(table(main.cluster.epi))

write.table(main.cluster.epi.counts, file = "moca_cluster_counts.txt", sep = "\t", quote = F)
