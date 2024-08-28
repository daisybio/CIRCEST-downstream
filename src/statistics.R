library(shiny)
library(bslib)
library(plotly)
library(heatmaply)
library(shinycssloaders)

statisticsUI <- function(id) {
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      card(
        card_header("Differential expression options"),
        card_body(
          uiOutput(ns("column_select")),
          uiOutput(ns("control_select")),
          uiOutput(ns("treatment_select")),
          class = "overflow"
        ),
        class = "overflow"
      ),
      input_task_button(ns("run_test"), "Run test")
    ),
    mainPanel(
      card(
        card_header("Plot options"),
        card_body(
          sliderInput(ns("alpha"),
            "Alpha",
            min = 0,
            max = 1,
            value = 0.05
          ),
          downloadButton(ns("download_correlation"), "Download results")
        )
      )
    )
  )
}

statisticsServer <- function(id, expr_matrix, phenotype) {
  moduleServer(id, function(input, output, session) {

    output$column_select <- renderUI({
      ns <- NS(id)
      selectInput(
        ns("col"),
        "Column",
        choices = colnames(phenotype()),
        selected = colnames(phenotype())[1]
      )
    })

    possible_vals <- reactive({
      req(input$col)
      unique(phenotype()[[input$col]])
    })

    output$control_select <- renderUI({
      ns <- NS(id)

      selectizeInput(
        ns("control"),
        "Control",
        multiple = TRUE,
        choices = possible_vals()
      )
    })

    output$treatment_select <- renderUI({
      ns <- NS(id)

      selectizeInput(
        ns("treatment"),
        "Treatment",
        multiple = TRUE,
        choices = possible_vals()
      )
    })
  })
}
