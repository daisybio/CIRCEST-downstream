library(shiny)
library(umap)
library(plotly)
library(bslib)
library(shinycssloaders)

dimredUI <- function(id) {
  ns <- NS(id)

  tagList(
    card(
      card_header("Settings"),
      card_body(
        uiOutput(ns("select_coloring")),
        class = "overflow"
      ),
      class = "overflow"
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

dimredServer <- function(id, expr_matrix, phenotype) {
  moduleServer(id, function(input, output, session) {
    output$select_coloring <- renderUI({
      ns <- NS(id)

      columns <- colnames(phenotype())

      selectInput(
        ns("coloring"),
        "Color by",
        choices = columns
      )
    })

    log_expression <- reactive({
      log1p(expr_matrix())
    })

    pca3 <- reactive({
      print("Calculating PCA3")
      pca <- prcomp(t(log_expression()), rank. = 3)
      components <- pca[["x"]]
      components <- data.frame(components)
      cbind(components, phenotype())
    })

    pca10 <- reactive({
      print("Calculating PCA10")
      prcomp(t(log_expression()), rank. = 10)
    })

    umap_data <- reactive({
      print("Calculating UMAP")
      data <- pca10()$x
      res <- umap(data,
        n_components = 3,
        n_neighbors = min(15, nrow(data)) - 1
      )
      layout <- res[["layout"]]
      layout <- data.frame(layout)
      cbind(layout, phenotype())
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
