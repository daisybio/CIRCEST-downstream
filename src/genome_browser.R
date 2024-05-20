library(shiny)
library(htmltools)

genomeBrowserUI <- function(id) {
  ns <- NS(id)

  uiOutput(ns("test"))
}

genomeBrowserServer <- function(id, genome) {
  moduleServer(id, function(input, output, session) {
    output$test <- renderUI({
      genome_string <- if (is.null(genome)) {
        ""
      } else {
        paste0("genome=", genome)
      }
      HTML(paste0("<genome-browser ", genome_string, "></genome-browser>"))
    })
  })
}
