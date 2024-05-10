library(shiny)
library(umap)
library(plotly)
library(bslib)
library(SummarizedExperiment)
library(shinycssloaders)

dimredUI <- function(id) {
  ns <- NS(id)

  tagList(
    card(
      card_header("Settings"),
      card_body(
        selectInput(
          ns("coloring"),
          "Color by:",
          choices = NULL
        )
      )
    ),
    card(
      card_header("PCA"),
      card_body(
        withSpinner(
          plotlyOutput(ns("plotPCA"))
        )
      )
    ),
    card(
      card_header("UMAP"),
      card_body(
        withSpinner(
          plotlyOutput(ns("plotUMAP"))
        )
      )
    )
  )
}

dimredServer <- function(id, filtered) {
  moduleServer(id, function(input, output, session) {
    pca3 <- reactive({
      print("Calculating PCA3")
      se <- filtered()
      pca <- prcomp(t(assay(se, "norm")), rank. = 3)
      components <- pca[["x"]]
      components <- data.frame(components)
      cbind(components, colData(filtered()))
    })

    pca10 <- reactive({
      print("Calculating PCA10")
      se <- filtered()
      prcomp(t(assay(se, "norm")), rank. = 10)
    })

    umap_data <- reactive({
      print("Calculating UMAP")
      data <- pca10()$x
      se.umap <- umap(data,
        n_components = 3,
        n_neighbors = min(15, nrow(data)) - 1
      )
      layout <- se.umap[["layout"]]
      layout <- data.frame(layout)
      cbind(layout, colData(filtered()))
    })

    observeEvent(filtered(), {
      updateSelectInput(session,
        "coloring",
        choices = colnames(colData(filtered()))
      )
    })

    output$plotPCA <- renderPlotly({
      print("Rendering PCA")
      coloring <- req(input$coloring)
      data <- pca3()

      if (!coloring %in% colnames(data)) {
        coloring <- colnames(data)[1]
      }

      plot_ly(
        data,
        x = ~PC1,
        y = ~PC2,
        z = ~PC3,
        color = ~ get(coloring),
        text = ~ paste("Sample: ", rownames(data))
      ) %>% add_markers()
    })

    output$plotUMAP <- renderPlotly({
      print("Rendering UMAP")
      coloring <- req(input$coloring)
      data <- umap_data()

      if (!coloring %in% colnames(data)) {
        coloring <- colnames(data)[1]
      }

      plot_ly(
        data,
        x = ~X1,
        y = ~X2,
        z = ~X3,
        color = ~ get(coloring),
        text = ~ paste("Sample: ", rownames(data))
      ) %>% add_markers()
    })
  })
}
