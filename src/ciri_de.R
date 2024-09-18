library(stringi)

work_dir <- "../data/work"
circ_dir <- "../data/circ"
stringtie_dir <- "../data/transcripts"

create_df <- function(samples, condition) {
  ciri_files <- file.path(circ_dir, paste0(samples, ".gtf"))
  stringtie_files <- file.path(stringtie_dir, paste0(samples, "_out.gtf"))
  conditions <- rep(condition, length(samples))
  return(data.frame("ciri" = ciri_files, "condition" = conditions, "stringtie" = stringtie_files, row.names = samples))
}

ciri_prepde <- function(df, wd) {
  samplesheet <- file.path(wd, "samplesheet_ciri.tsv")
  write.table(df, samplesheet, sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)

  lib_file <- file.path(wd, "lib_ciri.csv")
  circ_file <- file.path(wd, "circ_ciri.csv")
  bsj_file <- file.path(wd, "bsj_ciri.csv")
  ratio_file <- file.path(wd, "ratio_ciri.csv")

  system(paste("prep_CIRIquant -i", samplesheet, "--lib", lib_file, "--circ", circ_file, "--bsj", bsj_file, "--ratio", ratio_file))

  return(list(lib_file, circ_file, bsj_file, ratio_file))
}

stringtie_prepde <- function(df, wd) {
  samplesheet <- file.path(wd, "samplesheet_stringtie.tsv")
  write.table(df, samplesheet, sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)

  gene_file <- file.path(wd, "gene.csv")
  transcript_file <- file.path(wd, "transcript.csv")

  system(paste("prepDE.py -i", samplesheet, "-g", gene_file, "-t", transcript_file))

  return(list(gene_file, transcript_file))
}

ciriquant_de <- function(lib_file, bsj_file, gene_file, wd) {
  gene_results <- file.path(wd, "gene_results.csv")
  circ_results <- file.path(wd, "circ_results.csv")

  system(paste("CIRI_DE_replicate --lib", lib_file, "--bsj", bsj_file, "--gene", gene_file, "--out", circ_results, "--out2", gene_results))

  return(list(gene_results, circ_results))
}

run <- function(control, treatment) {
  control_df <- create_df(control, "C")
  treatment_df <- create_df(treatment, "T")
  df <- rbind(control_df, treatment_df)

  tryCatch({
    temp_dir <- file.path(work_dir, stri_rand_strings(1, 10))
    dir.create(temp_dir)

    p_ciri <- ciri_prepde(df[c("ciri", "condition")], temp_dir)
    p_stringtie <- stringtie_prepde(df[c("stringtie")], temp_dir)

    p_ciri_de <- ciriquant_de(p_ciri[[1]], p_ciri[[3]], p_stringtie[[1]], temp_dir)

    df_gene <- read.csv(p_ciri_de[[1]], row.names = 1)
    df_circ <- read.csv(p_ciri_de[[2]], row.names = 1)

    return(list(df_gene, df_circ))
  },
  finally = {
    unlink(temp_dir, recursive = TRUE)
  })
}
