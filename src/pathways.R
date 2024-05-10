library(shiny)
library(httr)
library(heatmaply)
library(plotly)
library(SummarizedExperiment)

organisms <- httr::content(
  GET("https://webservice.wikipathways.org/listOrganisms?format=json")
)$organisms

pathwaysUI <- function(id) {
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      selectInput(ns("organism"),
        "Organism",
        choices = organisms,
        selected = "Mus musculus"
      ),
      selectInput(ns("pathway"),
        "Pathway",
        choices = NULL
      ),
      uiOutput(ns("pathway_genes"))
    ),
    mainPanel(
      textOutput(ns("pathway_heatmap_alt")),
      plotlyOutput(ns("pathway_heatmap")),
      downloadButton(ns("download_pathway"), "Download SummarizedExperiment")
    )
  )
}

pathwaysServer <- function(id, filtered) {
  moduleServer(id, function(input, output, session) {
    pathways <- reactive({
      print("Getting pathways")
      httr::content(
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
      p <- data.frame(httr::content(GET(url)))

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

    observeEvent(pathway_names(), {
      print("Updating pathway names")
      updateSelectInput(session,
        "pathway",
        choices = pathway_names()
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
  })
}
