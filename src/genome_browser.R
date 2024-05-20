library(shiny)
library(htmltools)
library(SummarizedExperiment)
library(GenomicRanges)

genomeBrowserUI <- function(id) {
  ns <- NS(id)

  uiOutput(ns("genome_browser"))
}

genomeBrowserServer <- function(id, genome, se) {
  moduleServer(id, function(input, output, session) {
    output$genome_browser <- renderUI({
      genome_string <- if (is.null(genome)) {
        ""
      } else {
        paste0("genome=", genome)
      }
      granges <- rowData(se())$X
      locstrings <- paste0(
        seqnames(granges), ":",
        start(granges), "-",
        end(granges)
      )
      regions_string <- paste0("regions=", paste(locstrings, collapse = ","))
      HTML(paste0(
        "<genome-browser ",
        genome_string, " ", regions_string, "></genome-browser>"
      ))
    })
  })
}
