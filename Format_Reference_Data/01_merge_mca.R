library(Matrix)
setwd("~/Yan_SSD/R_Analysis/teratoma-analysis/Reference_Data/")

files <- list.files("rmbatch_dge/")
files <- paste0("rmbatch_dge/", files)

mca.counts <- as.matrix(read.table(files[[1]], sep = " ", header = T, row.names = 1))
mca.counts <- Matrix(mca.counts)

for (i in 2:length(files)) {
  print(i)
  fi.counts <- as.matrix(read.table(files[[i]], sep = " ", header = T, row.names = 1))
  fi.counts <- Matrix(fi.counts)
  
  merged.genes <- union(rownames(fi.counts), rownames(mca.counts))
  fi.pad.genes <- merged.genes[!merged.genes %in% rownames(fi.counts)]
  fi.pad.counts <- Matrix(0, nrow = length(fi.pad.genes), ncol = ncol(fi.counts),
                          dimnames = list(fi.pad.genes, colnames(fi.counts)))
  fi.counts <- rbind(fi.counts, fi.pad.counts)
  
  mca.pad.genes <- merged.genes[!merged.genes %in% rownames(mca.counts)]
  mca.pad.counts <- Matrix(0, nrow = length(mca.pad.genes), ncol = ncol(mca.counts),
                          dimnames = list(mca.pad.genes, colnames(mca.counts)))
  mca.counts <- rbind(mca.counts, mca.pad.counts)
  
  mca.counts <- cbind(mca.counts[merged.genes,], fi.counts[merged.genes,])
  print(dim(mca.counts))
}

print(dim(mca.counts))
save(mca.counts, file = "mca_full_debatched.RData")
