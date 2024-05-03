library(shiny)
library(bslib)
library(umap)
library(DESeq2)
library(ggfortify)
library(plotly)
library(fishpond)
library(heatmaply)

source("preprocessing.R")
source("statistics.R")
source("about.R")

# Define UI for app that draws a histogram ----
ui <- navbarPage(
  "circRNA investigator",
  preprocessingUI,
  statisticsUI,
  aboutUI,
  theme = bs_theme(version = 5, bootswatch = "shiny")
)

deseq_enabled <- FALSE

# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  dataset <- reactive({
    se <- readRDS(paste0(input$dataset, "/tx.rds"))
    rownames(colData(se)) <- colData(se)$names
    colData(se)$names <- NULL
    se
  })

  deseq_coldata <- reactive({
    se <- dataset()
    colData(se)[, sapply(
      colData(se),
      function(x) length(unique(x)) > 1
    )]
  })

  deseq_design <- reactive({
    colDataFiltered <- deseq_coldata()
    names <- colnames(colDataFiltered)
    formula(paste0("~", paste(names, collapse = "+")))
  })

  normalized_genes <- reactive({
    table <- read.table(paste0(input$dataset, "/gene.tsv"),
      header = TRUE, sep = "\t"
    )
    rownames(table) <- table$gene_name
    table$gene_name <- NULL
    table$gene_id <- NULL

    if (deseq_enabled) {
      dds <- DESeqDataSetFromMatrix(
        round(table, 0),
        deseq_coldata(),
        design = deseq_design()
      )
      dds <- DESeq(dds)
      normalized <- counts(dds, normalized = TRUE)
      log1p(normalized)
    } else {
      log1p(table)
    }
  })

  normalized <- reactive({
    se <- dataset()

    if (deseq_enabled) {
      # Normalize
      dds <- DESeqDataSetFromMatrix(
        round(assay(se), 0),
        deseq_coldata(),
        design = deseq_design()
      )
      dds <- DESeq(dds)

      assay(se, "norm") <- counts(dds, normalized = TRUE)
    } else {
      assay(se, "norm") <- assay(se)
    }

    se
  })

  filtered <- reactive({
    se <- normalized()

    n_samples <- ncol(se)
    se <- scaleInfReps(se)
    se <- labelKeep(
      se,
      input$min_count,
      n_samples * input$min_samples_pct / 100
    )
    se <- se[rowData(se)$keep, ]
    se
  })

  pca3 <- reactive({
    se <- filtered()
    pca <- prcomp(t(assay(se, "norm")), rank. = 3)
    components <- pca[["x"]]
    components <- data.frame(components)
    cbind(components, colData(filtered()))
  })

  pca10 <- reactive({
    se <- filtered()
    prcomp(t(assay(se, "norm")), rank. = 10)
  })

  umap_data <- reactive({
    data <- pca10()$x
    se.umap <- umap(data,
      n_components = 3,
      n_neighbors = min(15, nrow(data)) - 1
    )
    layout <- se.umap[["layout"]]
    layout <- data.frame(layout)
    cbind(layout, colData(filtered()))
  })

  output$colorings <- renderUI({
    selectInput("coloring",
      "Color by:",
      choices = colnames(colData(filtered()))
    )
  })

  output$plotPCA <- renderPlotly({
    req(input$coloring)
    data <- pca3()
    plot_ly(
      data,
      x = ~PC1,
      y = ~PC2,
      z = ~PC3,
      color = ~ get(input$coloring),
      text = ~ paste("Sample: ", rownames(data))
    ) %>% add_markers()
  })

  output$plotUMAP <- renderPlotly({
    req(input$coloring)
    data <- umap_data()
    plot_ly(
      data,
      x = ~X1,
      y = ~X2,
      z = ~X3,
      color = ~ get(input$coloring),
      text = ~ paste("Sample: ", rownames(data))
    ) %>% add_markers()
  })

  se_cor <- reactive({
    se <- filtered()
    colData(se) <- cbind(colData(se), t(normalized_genes()))
    se
  })

  output$statisticsUI <- renderUI({
    if (input$test_type == "correlation") {
      correlationUI
    } else {
      differentialUI
    }
  })

  cor_choices <- reactive({
    data <- se_cor()
    # Select all numeric columns
    col_numeric <- colnames(colData(data))[sapply(colData(data), is.numeric)]
    # Select all non-singular columns
    col_numeric[sapply(
      colData(data)[col_numeric],
      function(x) length(unique(x)) > 1
    )]
  })

  observeEvent(list(cor_choices(), input$test_type), {
    updateSelectizeInput(session,
      "cor_x",
      choices = cor_choices(), server = TRUE
    )
  })

  observeEvent(input$test_type, {
    columns <- colnames(colData(filtered()))
    # Keep columns with more than 1 unique value
    columns <- columns[sapply(
      colData(filtered())[columns],
      function(x) length(unique(x)) > 1
    )]

    updateSelectInput(session,
      "diffex_col",
      choices = columns
    )
  })

  observeEvent(input$diffex_col, {
    updateSelectInput(session,
      "diffex_a",
      choices = unique(colData(filtered())[[input$diffex_col]])
    )
  })

  observeEvent(input$diffex_col, {
    choices <- unique(colData(filtered())[[input$diffex_col]])

    updateSelectInput(session,
      "diffex_b",
      choices = choices,
      selected = choices[2]
    )
  })


  output$download_correlation <- downloadHandler(
    filename = function() {
      paste("correlation-",
        input$cor_x, "-", input$cor_type, ".tsv",
        sep = ""
      )
    },
    content = function(file) {
      write.table(
        assay(se_cor(), "norm"),
        col.names = NA,
        file = file,
        sep = "\t"
      )
    },
    contentType = "text/tab-separated-values"
  )

  cor_result <- eventReactive(input$run_correlation, {
    swish(se_cor(), input$cor_x, cor = input$cor_type)
  })

  output$cor_volcano <- renderPlotly({
    req(cor_result())
    data <- data.frame(rowData(cor_result()))
    data$significant <- ifelse(data$qvalue < input$cor_alpha,
      "significant", "not significant"
    )
    plot_ly(
      data,
      x = ~log2FC,
      y = ~ -log10(qvalue),
      color = ~significant,
      text = ~ paste(
        "transcript: ", rownames(data),
        "<br>qvalue: ", qvalue, "<br>log2FC: ", log2FC
      )
    ) %>% add_markers()
  })

  output$cor_heatmap <- renderPlotly({
    req(cor_result())
    se <- cor_result()
    se <- se[rowData(se)$qvalue < input$cor_alpha, ]

    # Stop if there are less than 2 transcripts
    if (nrow(se) < 2) {
      return(NULL)
    }

    heatmaply(
      assay(se, "norm"),
      xlab = "Samples",
      ylab = "Transcripts",
    )
  })
}

shinyApp(ui = ui, server = server)
