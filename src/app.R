library(shiny)
library(bslib)
library(umap)
library(ggfortify)
library(plotly)
library(fishpond)
library(heatmaply)
library(httr)

source("load.R")

source("preprocessing.R")
source("pathways.R")
source("statistics.R")
source("about.R")

# Define UI for app that draws a histogram ----
ui <- navbarPage(
  "circRNA investigator",
  preprocessingUI,
  pathwaysUI,
  statisticsUI,
  aboutUI,
  theme = bs_theme(version = 5, bootswatch = "shiny")
)

se <- loadTx()
normalized_genes <- loadGenes()


server <- function(input, output, session) {
  filtered <- filteringServer("filtering", se)

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

  output$colorings <- renderUI({
    print("Rendering colorings")
    req(filtered())
    selectInput("coloring",
      "Color by:",
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

  se_cor <- reactive({
    print("Creating SummarizedExperiment for correlation")
    se <- filtered()
    print(dim(colData(se)))
    print(dim(normalized_genes))
    colData(se) <- cbind(colData(se), t(normalized_genes))
    se
  })

  output$statisticsUI <- renderUI({
    print("Rendering statistics UI")
    if (input$test_type == "correlation") {
      correlationUI
    } else {
      differentialUI
    }
  })

  cor_choices <- reactive({
    print("Calculating correlation choices")
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
    print("Updating correlation choices")
    updateSelectizeInput(session,
      "cor_x",
      choices = cor_choices(), server = TRUE
    )
  })

  observeEvent(input$test_type, {
    print("Updating correlation type")
    filtered <- req(filtered())
    columns <- colnames(colData(filtered))
    # Keep columns with more than 1 unique value
    columns <- columns[sapply(
      colData(filtered)[columns],
      function(x) length(unique(x)) > 1
    )]

    updateSelectInput(session,
      "diffex_col",
      choices = columns
    )
  })

  observeEvent(input$diffex_col, {
    print("Updating differential expression groups A")
    updateSelectInput(session,
      "diffex_a",
      choices = unique(colData(filtered())[[input$diffex_col]])
    )
  })

  observeEvent(input$diffex_col, {
    print("Updating differential expression groups B")
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
        rowData(test_result()),
        col.names = NA,
        file = file,
        sep = "\t"
      )
    },
    contentType = "text/tab-separated-values"
  )

  test_result <- eventReactive(input$run_test, {
    print("Running test")
    if (input$test_type == "correlation") {
      swish(se_cor(), input$cor_x, cor = input$cor_type)
    } else {
      data <- filtered()
      data <- data[colData(data)[, input$diffex_col] %in% c(
        input$diffex_a, input$diffex_b
      ), ]
      swish(data, input$diffex_col)
    }
  })

  output$volcano <- renderPlotly({
    print("Rendering volcano plot")
    req(test_result())
    data <- data.frame(rowData(test_result()))
    data$significant <- ifelse(data$qvalue < input$alpha,
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

  output$heatmap <- renderPlotly({
    print("Rendering heatmap")
    req(test_result())
    se <- test_result()
    se <- se[rowData(se)$qvalue < input$alpha, ]

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

  pathways <- reactive({
    print("Getting pathways")
    content(
      GET("https://webservice.wikipathways.org/listPathways",
        query = list(
          format = "json",
          organism = input$organism
        )
      )
    )$pathways
  })

  pathway_names <- reactive({
    print("Getting pathway names")
    sapply(pathways(), function(x) x$name)
  })

  selected_pathway <- reactive({
    print("Getting selected pathway")
    available_pathways <- req(pathways())
    selected <- available_pathways[req(pathway_names()) == req(input$pathway)]

    selected <- selected[[1]]

    selected
  })

  participants <- reactive({
    print("Getting participants")
    selected <- req(selected_pathway())
    url <- paste0(
      "https://www.wikipathways.org/wikipathways-assets/pathways/",
      selected$id,
      "/",
      selected$id,
      "-datanodes.tsv"
    )
    p <- data.frame(content(GET(url)))

    unique(p[p$Type == "GeneProduct", ]$Label)
  })

  output$pathway_genes <- renderUI({
    print("Rendering pathway genes")
    # A button for each gene
    genes <- participants()
    lapply(genes, function(gene) {
      actionButton(
        gene,
        gene
      )
    })
  })

  output$select_pathway <- renderUI({
    print("Rendering select pathway")
    p_names <- pathway_names()

    if (length(p_names) == 0) {
      return(NULL)
    }

    selectInput("pathway",
      "Pathway",
      choices = p_names
    )
  })

  se_pathway <- reactive({
    print("Getting SummarizedExperiment for pathway")
    se <- filtered()
    genes <- participants()

    # If genes is NULL, return an empty SummarizedExperiment
    if (is.null(genes)) {
      return(SummarizedExperiment())
    }

    # rowData(se)$gene_name is a comma-separated list of gene names
    # We split it and check if any of the genes are in the pathway
    # Ignore case
    splitted <- strsplit(rowData(se)$gene_name, ",")
    keep <- sapply(splitted, function(x) any(tolower(x) %in% tolower(genes)))
    se[keep, ]
  })

  output$download_pathway <- downloadHandler(
    filename = function() {
      paste("pathway_", sub(" ", "-", req(input$pathway)), ".rds", sep = "")
    },
    content = function(file) {
      saveRDS(se_pathway(), file)
    },
  )

  output$pathway_heatmap <- renderPlotly({
    print("Rendering pathway heatmap")
    se <- se_pathway()

    if (nrow(se) == 0) {
      return(NULL)
    }

    heatmaply(
      assay(se, "norm"),
      xlab = "Samples",
      ylab = "Transcripts",
    )
  })

  output$pathway_heatmap_alt <- renderText({
    print("Rendering pathway heatmap alt text")
    se <- se_pathway()

    if (nrow(se) == 0) {
      return("No matching transcripts found")
    }

    paste(
      "Found expression information for", nrow(se),
      " transcripts from genes in the pathway"
    )
  })
}

shinyApp(ui = ui, server = server)
