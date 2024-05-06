library(shiny)
library(bslib)
library(umap)
library(DESeq2)
library(ggfortify)
library(plotly)
library(fishpond)
library(heatmaply)
library(httr)

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

deseq_enabled <- TRUE

# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  dataset <- reactive({
    se <- readRDS(paste0(input$dataset, "/tx.rds"))
    rownames(colData(se)) <- colData(se)$names
    colData(se)$names <- NULL
    rowData(se)$type <- ifelse(
      grepl("^circ_", rownames(se)),
      "circular",
      "linear"
    )
    se
  })

  deseq_coldata <- reactive({
    se <- dataset()
    colData(se)[sapply(
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
    rownames(table) <- table$gene_id
    table$gene_name <- NULL
    table$gene_id <- NULL

    if (deseq_enabled) {
      rounded <- round(table, 0)
      # Reorder columns to match colData
      rounded <- rounded[, rownames(deseq_coldata())]
      dds <- DESeqDataSetFromMatrix(
        rounded,
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

      assay(se, "norm") <- log1p(counts(dds, normalized = TRUE))
    } else {
      assay(se, "norm") <- log1p(assay(se))
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

    if (input$transcript_types != "both") {
      se <- se[rowData(se)$type == input$transcript_types, ]
    }

    se
  })

  output$filtered_description <- renderText({
    se <- filtered()
    paste(
      "Found", nrow(se), "matching transcripts"
    )
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

  test_result <- eventReactive(input$run_test, {
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
    sapply(pathways(), function(x) x$name)
  })

  selected_pathway <- reactive({
    available_pathways <- pathways()
    selected <- available_pathways[pathway_names() == input$pathway]

    if (length(selected) == 0) {
      return(NULL)
    }
    selected <- selected[[1]]

    selected
  })

  participants <- reactive({
    selected <- selected_pathway()
    if (is.null(selected)) {
      return(NULL)
    }
    url <- paste0(
      "https://www.wikipathways.org/wikipathways-assets/pathways/",
      selected$id,
      "/",
      selected$id,
      "-datanodes.tsv"
    )
    p <- data.frame(content(GET(url)))

    p[p$Type == "GeneProduct", ]$Label
  })

  output$pathway_genes <- renderUI({
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

  output$pathway_heatmap <- renderPlotly({
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
    se <- se_pathway()

    if (nrow(se) == 0) {
      return("No matching transcripts found")
    }

    paste("Found expression information for", nrow(se), " transcripts from genes in the pathway")
  })
}

shinyApp(ui = ui, server = server)
