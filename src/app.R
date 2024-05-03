library(shiny)
library(bslib)
library(umap)
library(DESeq2)
library(fishpond)
library(ggfortify)
library(plotly)

source("preprocessing.R")
source("about.R")

# Define UI for app that draws a histogram ----
ui <- navbarPage(
  "circRNA investigator",
  preprocessingUI,
  aboutUI
)

deseq_enabled <- FALSE

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  dataset <- reactive({
    se <- readRDS(paste0(input$dataset, "/tx.rds"))
    rownames(colData(se)) <- colData(se)$names
    colData(se)$names <- NULL
    se
  })

  normalized <- reactive({
    se <- dataset()

    if (deseq_enabled) {
      # Keep only columns with more than on unique value
      colDataFiltered <- colData(se)[, sapply(colData(se), function(x) length(unique(x)) > 1)]
      names <- colnames(colDataFiltered)
      design <- formula(paste0("~", paste(names, collapse = "+")))

      # Normalize
      dds <- DESeqDataSetFromMatrix(
        round(assay(se), 0),
        colDataFiltered,
        design = design
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

  pca2 <- reactive({
    se <- filtered()
    pca <- prcomp(t(assay(se, "norm")), rank. = 2)
    components <- pca[['x']]
    components <- data.frame(components)
    cbind(components, colData(filtered()))
  })

  pca10 <- reactive({
    se <- filtered()
    prcomp(t(assay(se, "norm")), rank. = 10)
  })

  umap_data <- reactive({
    se.umap <- umap(pca10()$x)
    layout <- se.umap[['layout']]
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
    plot_ly(
      pca2(),
      x = ~PC1,
      y = ~PC2,
      color = ~get(input$coloring),
      type = "scatter",
      mode = "markers"
    )
  })

  output$plotUMAP <- renderPlotly({
    plot_ly(
      umap_data(),
      x = ~X1,
      y = ~X2,
      color = ~get(input$coloring),
      type = "scatter",
      mode = "markers"
    )
  })

  output$datadescription <- renderPrint({
    filtered()
  })
}

shinyApp(ui = ui, server = server)
