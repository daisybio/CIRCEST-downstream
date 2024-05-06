library(SummarizedExperiment)
library(umap)

se <- readRDS("data/test/tx.rds")

print(rowData(se))
