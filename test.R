source("load.R")
source("ciri_de.R")


phenotype <- loadPhenotype()

control_samples <- rownames(head(phenotype))
treatment_samples <- rownames(tail(phenotype))

results <- run(control_samples, treatment_samples)

print(head(results[[2]]))
