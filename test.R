library(SummarizedExperiment)
library(umap)

se <- readRDS("data/antiHormonal/tx.rds")

umap(t(assay(se, "counts")))

print("Finished")