library(shiny)
library(bslib)
library(ggfortify)
library(plotly)
library(fishpond)
library(heatmaply)
library(httr)

source("load.R")

source("filtering.R")
source("dimred.R")
source("pathways.R")
source("statistics.R")
source("about.R")

# Define UI for app that draws a histogram ----
ui <- navbarPage(
  "circRNA investigator",
  tabPanel(
    "Preprocessing",
    sidebarLayout(
      sidebarPanel(
        filteringUI("filtering")
      ),
      mainPanel(
        dimredUI("dimred")
      )
    )
  ),
  tabPanel(
    "Pathways",
    pathwaysUI("pathways"),
  ),
  tabPanel(
    "Statistics",
    statisticsUI("statistics"),
  ),
  aboutUI,
  theme = bs_theme(version = 5, bootswatch = "shiny")
)

se <- loadTx()
normalized_genes <- loadGenes()


server <- function(input, output, session) {
  filtered <- filteringServer("filtering", se)
  dimredServer("dimred", filtered)
  pathwaysServer("pathways", filtered)
  statisticsServer("statistics", filtered, normalized_genes)
}

shinyApp(ui = ui, server = server)
