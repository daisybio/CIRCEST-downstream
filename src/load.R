data_prefix <- "../data/"

loadGenes <- function() {
  genes_path <- paste0(data_prefix, "gene.joined.tsv")

  genes <- read.table(genes_path, header = TRUE, sep = "\t")
  rownames(genes) <- genes$gene_id
  genes$gene_id <- NULL
  genes$gene_name <- NULL

  return(genes)
}

loadCirc <- function() {
  circ_path <- paste0(data_prefix, "circ.joined.tsv")

  circ <- read.table(circ_path, header = TRUE, sep = "\t")
  rownames(circ) <- circ$circ_id
  circ$circ_id <- NULL
  circ$gene_name <- NULL

  return(circ)
}

loadPhenotype <- function() {
  phenotype_path <- paste0(data_prefix, "phenotype.csv")

  phenotype <- read.csv(phenotype_path, header = TRUE, sep = ",")
  rownames(phenotype) <- phenotype$sample
  phenotype$sample <- NULL

  return(phenotype)
}

loadGenome <- function() {
  # Return the content of genome.txt if it exists
  genome_path <- paste0(data_prefix, "genome.txt")
  if (file.exists(genome_path)) {
    return(readLines(genome_path))
  }
}
