library(shiny)
library(bslib)

source("load.R")

source("filter_transcripts.R")
source("dimred.R")
source("pathways.R")
source("statistics.R")
source("about.R")

# Define UI for app that draws a histogram ----
ui <- navbarPage(
  "circRNA investigator",
  tags$head(tags$style(HTML("
    .overflow {
      overflow: visible !important;
    }
  "))),
  tabPanel(
    "Preprocessing",
    sidebarLayout(
      sidebarPanel(
        filterTranscriptsUI("filtering")
      ),
      mainPanel(
        dimredUI("dimred")
      )
    )
  ),
  tabPanel(
    "Pathways",
    pathwaysUI("pathways")
  ),
  tabPanel(
    "Statistics",
    statisticsUI("statistics")
  ),
  tabPanel(
    "About",
    aboutUI
  ),
  id = "navbar",
  theme = bs_theme(version = 5, bootswatch = "shiny")
)

se <- loadTx()
normalized_genes <- loadGenes()

server <- function(input, output, session) {
  filtered <- filterTranscriptsServer("filtering", se)
  dimredServer("dimred", filtered, colnames(colData(se)))
  pathwaysServer("pathways", filtered)
  statisticsServer(
    "statistics", filtered, normalized_genes,
    reactive(input$navbar)
  )
}

shinyApp(ui = ui, server = server)
