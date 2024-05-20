library(shiny)
library(htmltools)

genomeBrowserUI <- function(id) {
  ns <- NS(id)

  uiOutput(ns("test"))
}

genomeBrowserServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    output$test <- renderUI({
      HTML("<genome-browser></genome-browser>")
    })
  })
}
