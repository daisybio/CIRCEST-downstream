library(shiny)
library(httr)

organisms <- content(
  GET("https://webservice.wikipathways.org/listOrganisms?format=json")
)$organisms

pathwaysUI <- tabPanel(
  "Pathways",
  sidebarLayout(
    sidebarPanel(
      selectInput("organism",
        "Organism",
        choices = organisms,
        selected = "Mus musculus"
      ),
      uiOutput("select_pathway"),
      uiOutput("pathway_genes")
    ),
    mainPanel(
      textOutput("pathway_heatmap_alt"),
      plotlyOutput("pathway_heatmap")
    )
  )
)
