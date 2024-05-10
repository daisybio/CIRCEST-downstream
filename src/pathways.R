library(shiny)
library(httr)
library(heatmaply)
library(plotly)
library(SummarizedExperiment)
library(bslib)
library(shinycssloaders)

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
      withSpinner(
        uiOutput(ns("pathway_selector")),
        proxy.height = "50px"
      ),
      withSpinner(
        uiOutput(ns("pathway_genes"))
      )
    ),
    mainPanel(
      card(
        card_header("Downloads"),
        card_body(
          uiOutput(ns("download_buttons"))
        )
      ),
      card(
        card_header("Heatmap"),
        card_body(
          textOutput(ns("pathway_heatmap_alt")),
          withSpinner(
            plotlyOutput(ns("pathway_heatmap"))
          )
        )
      )
    )
  )
}

pathwaysServer <- function(id, filtered) {
  moduleServer(id, function(input, output, session) {
    pathways <- reactive({
      print("Getting pathways")
    })

    output$pathway_selector <- renderUI({
      pathways <- httr::content(
        GET("https://webservice.wikipathways.org/listPathways",
          query = list(
            format = "json",
            organism = input$organism
          )
        )
      )$pathways

      ns <- NS(id)

      pathwayIDs <- sapply(pathways, function(x) x$id)
      names(pathwayIDs) <- sapply(pathways, function(x) x$name)

      selectInput(
        ns("pathway"),
        "Pathway",
        choices = pathwayIDs
      )
    })

    participants <- reactive({
      print("Getting participants")
      req(input$organism)
      selected <- req(input$pathway)
      url <- paste0(
        "https://www.wikipathways.org/wikipathways-assets/pathways/",
        selected,
        "/",
        selected,
        "-datanodes.tsv"
      )
      p <- data.frame(httr::content(GET(url)))

      unique(p[p$Type == "GeneProduct", ]$Label)
    })

    output$pathway_genes <- renderUI({
      print("Rendering pathway genes")
      # A button for each gene
      genes <- participants()
      div(
        lapply(genes, function(gene) {
          actionButton(
            gene,
            gene,
            style = "flex-grow: 1; margin: 1px"
          )
        }),
        style = "display: flex; flex-wrap: wrap; "
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

    output$download_buttons <- renderUI({
      ns <- NS(id)
      print("Rendering download buttons")
      req(input$organism)
      req(input$pathway)
      req(se_pathway())

      div(
        downloadButton(
          ns("download_csv"),
          "Download CSV"
        ),
        downloadButton(
          ns("download_se"),
          "Download SummarizedExperiment"
        )
      )
    })

    output$download_csv <- downloadHandler(
      filename = function() {
        paste("pathway_", sub(" ", "-", req(input$pathway)), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(assay(se_pathway()), file, row.names = TRUE)
      },
    )

    output$download_se <- downloadHandler(
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
