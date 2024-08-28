source("load.R")
library(umap)


circ <- loadCirc()
phenotype <- loadPhenotype()

circ <- log1p(circ)

pca10 <- prcomp(t(circ), rank. = 10)

data <- pca10$x

res <- umap(data,
  n_components = 3,
  n_neighbors = min(15, nrow(data)) - 1
)

layout <- res[["layout"]]
layout <- data.frame(layout)

print(dim(layout))

cbind(layout, phenotype)