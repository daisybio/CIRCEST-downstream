library(SummarizedExperiment)

data_prefix <- "../data/"

se <- readRDS(paste0(data_prefix, "tx.rds"))
rownames(colData(se)) <- colData(se)$names
colData(se)$names <- NULL
rowData(se)$type <- ifelse(
  grepl("^circ_", rownames(se)),
  "circular",
  "linear"
)
assay(se, "log") <- log1p(assay(se, "tpm"))

genes_table <- read.table(paste0(data_prefix, "gene.tsv"),
  header = TRUE, sep = "\t"
)
# Remove leading X characters from colnames
colnames(genes_table) <- gsub("^X", "", colnames(genes_table))
rownames(genes_table) <- genes_table$gene_id
genes_table$gene_name <- NULL
genes_table$gene_id <- NULL
genes_log <- log1p(genes_table)

loadTx <- function() {
  return(se)
}

loadGenes <- function() {
  return(genes_log)
}