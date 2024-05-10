library(shiny)
library(bslib)

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
  theme = bs_theme(version = 5, bootswatch = "shiny")
)

se <- loadTx()
normalized_genes <- loadGenes()

server <- function(input, output, session) {
  filtered <- filteringServer("filtering", se)
  dimredServer("dimred", filtered, colnames(colData(se)))
  pathwaysServer("pathways", filtered)
  statisticsServer("statistics", filtered, normalized_genes)
}

shinyApp(ui = ui, server = server)
