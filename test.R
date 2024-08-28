source("load.R")

phenotype <- loadPhenotype()

keep_values <- unique(phenotype[['ag']])
keep_values <- as.character(keep_values)

print(phenotype[['age']] %in% keep_values)
