library(shiny)
library(bslib)
library(shinycssloaders)

filterTranscriptsUI <- function(id) {
  ns <- NS(id)

  tagList(
    card(
      card_header("Filter transcripts"),
      card_body(
        uiOutput(ns("threshold_selector")),
        sliderInput(ns("min_samples_pct"),
          "Min samples %",
          min = 0,
          max = 100,
          value = 20.0
        ),
        textOutput(ns("filtered_description"))
      )
    ),
    withSpinner(
      uiOutput(ns("download_button")),
      proxy.height = "50px"
    )
  )
}

filterTranscriptsServer <- function(id, matrix, phenotype) {
  moduleServer(id, function(input, output, session) {
    matrixSamples <- reactive({
      cur_matrix <- matrix[, rownames(phenotype())]
      cur_matrix$gene_id <- NULL
      cur_matrix
    })

    output$threshold_selector <- renderUI({
      cur_matrix <- matrixSamples()
      ns <- NS(id)

      sliderInput(ns("threshold"),
        "Threshold",
        min = 0,
        max = max(cur_matrix),
        value = 0
      )
    })

    filtered <- reactive({
      cur_matrix <- matrixSamples()
      threshold <- input$threshold

      keep <- cur_matrix >= threshold
      keep <- rowSums(keep) >= input$min_samples_pct / 100 * ncol(cur_matrix)

      cur_matrix[keep, ]
    })

    output$filtered_description <- renderText({
      paste(
        "Found", nrow(filtered()), "matching transcripts"
      )
    })

    output$download_button <- renderUI({
      ns <- NS(id)
      req(filtered())
      downloadButton(
        ns("download_filtered"),
        "Download filtered expression matrix"
      )
    })

    output$download_filtered <- downloadHandler(
      filename = function() {
        "filtered_expression.tsv"
      },
      content = function(file) {
        write.table(filtered(), file, sep = "\t", col.names = NA)
      }
    )

    return(filtered)
  })
}
