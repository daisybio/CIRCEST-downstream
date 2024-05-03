library(shiny)
library(bslib)
library(umap)
library(DESeq2)
library(ggfortify)
library(plotly)
library(fishpond)

source("preprocessing.R")
source("correlation.R")
source("about.R")

# Define UI for app that draws a histogram ----
ui <- navbarPage(
  "circRNA investigator",
  preprocessingUI,
  correlationUI,
  aboutUI
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
      counts(dds, normalized = TRUE)
    } else {
      table
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

  observeEvent(cor_choices(), {
    updateSelectizeInput(session,
      "cor_x",
      choices = cor_choices(), server = TRUE
    )
  })

  cor_result <- eventReactive(input$run_correlation, {
    swish(se_cor(), input$cor_x, cor = input$cor_type)
  })

  output$cor_volcano <- renderPlotly({
    req(cor_result())
    data <- data.frame(rowData(cor_result()))
    plot_ly(
      data,
      x = ~log2FC,
      y = ~ -log10(qvalue),
      color = ~ ifelse(qvalue < 0.05, "red", "black"),
      text = ~ paste("Gene: ", rownames(data))
    ) %>% add_markers()
  })

  output$summary <- renderPrint({
    req(cor_result())
    rowData(cor_result())
  })
}

shinyApp(ui = ui, server = server)
