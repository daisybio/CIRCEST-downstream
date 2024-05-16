library(shiny)
library(fishpond)
library(SummarizedExperiment)
library(bslib)
library(shinycssloaders)

filterSamplesUI <- function(id) {
  ns <- NS(id)

  card(
    card_header("Hello"),
    card_body(
      "world"
    )
  )
}

filterSamplesServer <- function(id, se) {
  moduleServer(id, function(input, output, session) {
    se <- se

    filtered <- reactive({
      se
    })

    return(filtered)
  })
}