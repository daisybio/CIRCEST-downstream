library(shiny)
library(bslib)

source("load.R")

source("filter_samples.R")
source("filter_transcripts.R")
source("dimred.R")
source("pathways.R")
source("statistics.R")
source("genome_browser.R")
source("about.R")

# Define UI for app that draws a histogram ----
ui <- navbarPage(
  "circRNA investigator",
  tags$head(tags$style(HTML("
    .overflow {
      overflow: visible !important;
    }
  "))),
  tags$script(src = "bundle.js"),
  tabPanel(
    "Preprocessing",
    sidebarLayout(
      sidebarPanel(
        filterSamplesUI("filter_samples"),
        filterTranscriptsUI("filter_transcripts")
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
    "Genome browser",
    genomeBrowserUI("genome_browser")
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
genome <- loadGenome()

server <- function(input, output, session) {
  se_filtered_samples <- filterSamplesServer("filter_samples", se)
  filtered <- filterTranscriptsServer("filter_transcripts", se_filtered_samples)
  dimredServer("dimred", filtered, colnames(colData(se)))
  pathwaysServer("pathways", filtered)
  statisticsServer(
    "statistics", filtered, normalized_genes,
    reactive(input$navbar)
  )
  genomeBrowserServer("genome_browser", genome)
}

shinyApp(ui = ui, server = server)
