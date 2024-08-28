library(shiny)
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

filterSamplesServer <- function(id, phenotype) {
  moduleServer(id, function(input, output, session) {
    output$columnsUI <- renderUI({
      ns <- NS(id)
      columns <- colnames(phenotype)
      inputs <- lapply(columns, function(column) {
        selectizeInput(
          inputId = ns(column),
          label = column,
          choices = unique(phenotype[column]),
          multiple = TRUE
        )
      })
      tagList(inputs)
    })

    filtered <- reactive({
      print("Filtering samples")
      columns <- colnames(phenotype)

      keep <- rep(TRUE, nrow(phenotype))

      for (i in seq_along(columns)) {
        column <- columns[i]
        selected <- input[[column]]
        if (is.null(selected)) {
          selected <- unique(phenotype[[column]])
        }
        selected <- as.character(selected)
        current <- phenotype[[column]] %in% selected
        keep <- keep & current
      }
      subset(phenotype, keep)
    })

    output$filtered_description <- renderText({
      phenotype_filtered <- filtered()
      paste0("Found ", nrow(phenotype_filtered), " matching samples")
    })

    return(filtered)
  })
}
