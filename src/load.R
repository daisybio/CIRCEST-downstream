library(DESeq2)

deseq_enabled <- TRUE

data_prefix <- "../data/"

se <- readRDS(paste0(data_prefix, "tx.rds"))
rownames(colData(se)) <- colData(se)$names
colData(se)$names <- NULL
rowData(se)$type <- ifelse(
  grepl("^circ_", rownames(se)),
  "circular",
  "linear"
)

deseq_coldata <- colData(se)[sapply(
  colData(se),
  function(x) length(unique(x)) > 1
)]

deseq_design <- formula(paste0(
  "~",
  paste(colnames(deseq_coldata), collapse = "+")
))

genes_table <- read.table(paste0(data_prefix, "gene.tsv"),
  header = TRUE, sep = "\t"
)
# Remove leading X characters from colnames
colnames(genes_table) <- gsub("^X", "", colnames(genes_table))
rownames(genes_table) <- genes_table$gene_id
genes_table$gene_name <- NULL
genes_table$gene_id <- NULL

loadTx <- function() {
  if (deseq_enabled) {
    dds <- DESeqDataSetFromMatrix(
      round(assay(se), 0),
      deseq_coldata,
      design = deseq_design
    )

    dds <- DESeq(dds)

    assay(se, "norm") <- log1p(counts(dds, normalized = TRUE))
  } else {
    assay(se, "norm") <- log1p(assay(se))
  }

  return(se)
}

loadGenes <- function() {
  if (deseq_enabled) {
    rounded <- round(genes_table, 0)
    # Reorder columns to match colData
    rounded <- rounded[, rownames(deseq_coldata)]
    dds <- DESeqDataSetFromMatrix(
      rounded,
      deseq_coldata,
      design = deseq_design
    )
    dds <- DESeq(dds)
    normalized_genes <- counts(dds, normalized = TRUE)
    normalized_genes <- log1p(normalized_genes)
  } else {
    normalized_genes <- log1p(genes_table)
  }

  return(normalized_genes)
}