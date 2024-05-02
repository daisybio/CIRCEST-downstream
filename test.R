library(SummarizedExperiment)

se <- readRDS("data/antiHormonal/tx.rds")

rownames(colData(se)) <- colData(se)$names
colData(se)$names <- NULL

print(colData(se))
