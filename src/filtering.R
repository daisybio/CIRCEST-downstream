library(shiny)
library(fishpond)
library(SummarizedExperiment)
library(bslib)

filteringUI <- function(id) {
  ns <- NS(id)

  card(
    card_header("Filter transcripts"),
    card_body(
      sliderInput(ns("min_tpm"),
        "Min TPM",
        min = 0,
        max = 100,
        value = 20
      ),
      sliderInput(ns("min_samples_pct"),
        "Min samples %",
        min = 0,
        max = 100,
        value = 20.0
      ),
      radioButtons(ns("transcript_types"), "Transcript types",
        choices = c("circular", "linear", "both"),
        selected = "circular"
      ),
      textOutput(ns("filtered_description"))
    )
  )
}

labelKeepTPM <- function(y, minTPM = 10, minN = 3, x) {
  if (!missing(x)) {
    stopifnot(x %in% names(colData(y)))
    minN <- min(table(colData(y)[[x]]))
    if (minN > 10) {
      minN <- 10 + (minN - 10) * 0.7
    }
  }
  cts <- assays(y)[["tpm"]]
  if (is(cts, "dgCMatrix")) {
    keep <- Matrix::rowSums(cts >= minTPM) >= minN
  } else {
    keep <- rowSums(cts >= minTPM) >= minN
  }
  mcols(y)$keep <- keep
  metadata(y)$preprocessed <- TRUE
  if (!"infRepsScaled" %in% names(metadata(y))) {
    metadata(y)$infRepsScaled <- FALSE
  }
  y
}

filteringServer <- function(id, se) {
  moduleServer(id, function(input, output, session) {
    filtered <- reactive({
      n_samples <- ncol(se)
      se <- scaleInfReps(se)
      se <- labelKeepTPM(
        se,
        input$min_tpm,
        n_samples * input$min_samples_pct / 100
      )
      se <- se[rowData(se)$keep, ]

      if (input$transcript_types != "both") {
        se <- se[rowData(se)$type == input$transcript_types, ]
      }
      se
    })

    output$filtered_description <- renderText({
      paste(
        "Found", nrow(filtered()), "matching transcripts"
      )
    })

    return(filtered)
  })
}
