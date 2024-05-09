library(shiny)

filteringUI <- function(id) {
  ns <- NS(id)

  card(
    card_header("Filter transcripts"),
    card_body(
      sliderInput(ns("min_count"),
        "Min count",
        min = 0,
        max = 100,
        value = 50
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

filteringServer <- function(id, se) {
  moduleServer(id, function(input, output, session) {
    filtered <- reactive({
      n_samples <- ncol(se)
      se <- scaleInfReps(se)
      se <- labelKeep(
        se,
        input$min_count,
        n_samples * input$min_samples_pct / 100
      )
      se <- se[rowData(se)$keep, ]

      if (input$transcript_types != "both") {
        se <- se[rowData(se)$type == input$transcript_types, ]
      }
    })

    output$filtered_description <- renderText({
      paste(
        "Found", nrow(filtered()), "matching transcripts"
      )
    })

    return(filtered)
  })
}
