library(shiny)
library(bslib)
library(plotly)
library(heatmaply)
library(shinycssloaders)
source("ciri_de.R")

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
      ),
      card(
        card_header("Results: circRNA"),
        card_body(
          tableOutput(ns("result_table_circ"))
        )
      ),
      card(
        card_header("Results: gene"),
        card_body(
          tableOutput(ns("result_table_gene"))
        )
      )
    )
  )
}

statisticsServer <- function(id, expr_matrix, phenotype, circ_cpm) {
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

    test_result <- eventReactive(input$run_test, {
      print("Running test")
      column <- input$col
      control <- input$control
      treatment <- input$treatment

      cur_phenotype <- phenotype()

      control_samples <- rownames(cur_phenotype[cur_phenotype[[column]] %in% control, ])
      treatment_samples <- rownames(cur_phenotype[cur_phenotype[[column]] %in% treatment, ])

      run(control_samples, treatment_samples)
    })

    gene_results <- reactive({
      results <- test_result()[[1]]
      results[results$FDR < input$alpha, ]
    })

    circ_results <- reactive({
      results <- test_result()[[2]]
      results <- results[results$FDR < input$alpha, ]
      annotated <- cbind(circ_cpm[rownames(results), ]$gene_id, results)
      colnames(annotated)[1] <- "host_gene"
      annotated
    })

    output$result_table_circ <- renderTable({
      results <- circ_results()
      cbind("circRNA" = rownames(results), results)
    })

    output$result_table_gene <- renderTable({
      results <- gene_results()
      cbind("gene" = rownames(results), results)
    })
  })
}
