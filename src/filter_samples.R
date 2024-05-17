library(shiny)
library(fishpond)
library(SummarizedExperiment)
library(bslib)
library(shinycssloaders)

filterSamplesUI <- function(id) {
  ns <- NS(id)

  card(
    card_header("Filter samples"),
    card_body(
      uiOutput(ns("columnsUI")),
      textOutput(ns("filtered_description")),
      class = "overflow"
    ),
    class = "overflow"
  )
}

filterSamplesServer <- function(id, se) {
  moduleServer(id, function(input, output, session) {
    output$columnsUI <- renderUI({
      ns <- NS(id)
      columns <- colnames(colData(se))
      inputs <- lapply(columns, function(column) {
        selectizeInput(
          inputId = ns(column),
          label = column,
          choices = unique(colData(se)[[column]]),
          multiple = TRUE
        )
      })
      tagList(inputs)
    })

    filtered <- reactive({
      print("Filtering samples")
      columns <- colnames(colData(se))
      filters <- lapply(columns, function(column) {
        selected <- input[[column]]
        if (is.null(selected)) {
          selected <- unique(colData(se)[[column]])
        }
        selected
      })
      print(filters)

      keep <- rep(TRUE, nrow(colData(se)))
      for (i in seq_along(columns)) {
        keep <- keep & colData(se)[[columns[i]]] %in% filters[[i]]
      }
      se_filtered <- se[, keep]
      se_filtered
    })

    output$filtered_description <- renderText({
      se_filtered <- filtered()
      paste0("Found ", nrow(colData(se_filtered)), " matching samples")
    })

    return(filtered)
  })
}