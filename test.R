library(SummarizedExperiment)
library(umap)

se <- readRDS("data/antiHormonal/tx.rds")

print(colData(se))
