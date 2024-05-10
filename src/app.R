library(shiny)
library(bslib)
library(ggfortify)
library(plotly)
library(fishpond)
library(heatmaply)
library(httr)

source("load.R")

source("filtering.R")
source("dimred.R")
source("pathways.R")
# source("statistics.R")
source("about.R")

# Define UI for app that draws a histogram ----
ui <- navbarPage(
  "circRNA investigator",
  tabPanel(
    "Preprocessing",
    sidebarLayout(
      sidebarPanel(
        filteringUI("filtering")
      ),
      mainPanel(
        dimredUI("dimred")
      )
    )
  ),
  tabPanel(
    "Pathways",
    pathwaysUI("pathways"),
  ),
  # statisticsUI,
  aboutUI,
  theme = bs_theme(version = 5, bootswatch = "shiny")
)

se <- loadTx()
normalized_genes <- loadGenes()


server <- function(input, output, session) {
  filtered <- filteringServer("filtering", se)
  dimredServer("dimred", filtered)
  pathwaysServer("pathways", filtered)

  #  se_cor <- reactive({
  #    print("Creating SummarizedExperiment for correlation")
  #    se <- filtered()
  #    print(dim(colData(se)))
  #    print(dim(normalized_genes))
  #    colData(se) <- cbind(colData(se), t(normalized_genes))
  #    se
  #  })
  #
  #  output$statisticsUI <- renderUI({
  #    print("Rendering statistics UI")
  #    if (input$test_type == "correlation") {
  #      correlationUI
  #    } else {
  #      differentialUI
  #    }
  #  })
  #
  #  cor_choices <- reactive({
  #    print("Calculating correlation choices")
  #    data <- se_cor()
  #    # Select all numeric columns
  #    col_numeric <- colnames(colData(data))[sapply(colData(data), is.numeric)]
  #    # Select all non-singular columns
  #    col_numeric[sapply(
  #      colData(data)[col_numeric],
  #      function(x) length(unique(x)) > 1
  #    )]
  #  })
  #
  #  observeEvent(list(cor_choices(), input$test_type), {
  #    print("Updating correlation choices")
  #    updateSelectizeInput(session,
  #      "cor_x",
  #      choices = cor_choices(), server = TRUE
  #    )
  #  })
  #
  #  observeEvent(input$test_type, {
  #    print("Updating correlation type")
  #    filtered <- req(filtered())
  #    columns <- colnames(colData(filtered))
  #    # Keep columns with more than 1 unique value
  #    columns <- columns[sapply(
  #      colData(filtered)[columns],
  #      function(x) length(unique(x)) > 1
  #    )]
  #
  #    updateSelectInput(session,
  #      "diffex_col",
  #      choices = columns
  #    )
  #  })
  #
  #  observeEvent(input$diffex_col, {
  #    print("Updating differential expression groups A")
  #    updateSelectInput(session,
  #      "diffex_a",
  #      choices = unique(colData(filtered())[[input$diffex_col]])
  #    )
  #  })
  #
  #  observeEvent(input$diffex_col, {
  #    print("Updating differential expression groups B")
  #    choices <- unique(colData(filtered())[[input$diffex_col]])
  #
  #    updateSelectInput(session,
  #      "diffex_b",
  #      choices = choices,
  #      selected = choices[2]
  #    )
  #  })
  #
  #
  #  output$download_correlation <- downloadHandler(
  #    filename = function() {
  #      paste("correlation-",
  #        input$cor_x, "-", input$cor_type, ".tsv",
  #        sep = ""
  #      )
  #    },
  #    content = function(file) {
  #      write.table(
  #        rowData(test_result()),
  #        col.names = NA,
  #        file = file,
  #        sep = "\t"
  #      )
  #    },
  #    contentType = "text/tab-separated-values"
  #  )
  #
  #  test_result <- eventReactive(input$run_test, {
  #    print("Running test")
  #    if (input$test_type == "correlation") {
  #      swish(se_cor(), input$cor_x, cor = input$cor_type)
  #    } else {
  #      data <- filtered()
  #      data <- data[colData(data)[, input$diffex_col] %in% c(
  #        input$diffex_a, input$diffex_b
  #      ), ]
  #      swish(data, input$diffex_col)
  #    }
  #  })
  #
  #  output$volcano <- renderPlotly({
  #    print("Rendering volcano plot")
  #    req(test_result())
  #    data <- data.frame(rowData(test_result()))
  #    data$significant <- ifelse(data$qvalue < input$alpha,
  #      "significant", "not significant"
  #    )
  #    plot_ly(
  #      data,
  #      x = ~log2FC,
  #      y = ~ -log10(qvalue),
  #      color = ~significant,
  #      text = ~ paste(
  #        "transcript: ", rownames(data),
  #        "<br>qvalue: ", qvalue, "<br>log2FC: ", log2FC
  #      )
  #    ) %>% add_markers()
  #  })
  #
  #  output$heatmap <- renderPlotly({
  #    print("Rendering heatmap")
  #    req(test_result())
  #    se <- test_result()
  #    se <- se[rowData(se)$qvalue < input$alpha, ]
  #
  #    # Stop if there are less than 2 transcripts
  #    if (nrow(se) < 2) {
  #      return(NULL)
  #    }
  #
  #    heatmaply(
  #      assay(se, "norm"),
  #      xlab = "Samples",
  #      ylab = "Transcripts",
  #    )
  #  })
  #

  #
}

shinyApp(ui = ui, server = server)
