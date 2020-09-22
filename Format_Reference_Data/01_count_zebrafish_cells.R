cluster_id_map <- read.table("Reference_Data/Zebrafish/GSE112294_ClusterNames.txt", sep = "\t", header = T)
rownames(cluster_id_map) <- cluster_id_map$ClusterID

files.use <- list.files(path = "Reference_Data/Zebrafish/", pattern = "*.clustID.txt.gz")
clusters.list <- lapply(files.use, function(fi) {
  x <- scan(paste0("Reference_Data/Zebrafish/", fi), what = character())
  if(length(x) > 0) return(as.character(cluster_id_map[x, "TissueName"]))
  else return(c())
})


clusters <- unlist(clusters.list, F, F)
clusters[grepl("Neural", clusters)] <- "Ectoderm"
clusters[grepl("Epidermal", clusters)] <- "Ectoderm"

table(clusters)

