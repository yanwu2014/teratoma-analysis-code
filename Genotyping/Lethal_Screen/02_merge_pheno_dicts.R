library(methods)
library(perturbLM)

# Dataset names specified via command line
options <- commandArgs(trailingOnly = T)
output.name <- options[[1]]
input.files <- options[2:length(options)]

# Load genotype dictionaries
pheno.dicts <- lapply(1:length(input.files), function(i) {
  genotypes.list <- ReadGenotypes(input.files[[i]])
  lapply(genotypes.list, function(cells) paste(cells, i, sep = "_"))
})

genotypes <- names(pheno.dicts[[1]])
for (i in 2:length(pheno.dicts)) {
  genotypes <- intersect(genotypes, names(pheno.dicts[[i]]))
}
print(genotypes)

merged.dict <- vector("list", length(genotypes))
names(merged.dict) <- genotypes

for (dict in pheno.dicts) {
  for (g in genotypes) {
    merged.dict[[g]] <- c(merged.dict[[g]], dict[[g]])
  }
}

print(sapply(merged.dict, length))
print(sum(sapply(merged.dict, length)))

WriteGenotypes(merged.dict, output.name)
