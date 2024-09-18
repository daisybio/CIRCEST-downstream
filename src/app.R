library(shiny)
library(bslib)

source("load.R")

source("filter_samples.R")
source("filter_transcripts.R")
source("dimred.R")
source("statistics.R")
source("about.R")

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
    "Differential expression",
    statisticsUI("statistics")
  ),
  id = "navbar",
  theme = bs_theme(version = 5, bootswatch = "shiny")
)

phenotype <- loadPhenotype()
circ_cpm <- loadCirc()
gene_tpm <- loadGenes()
genome <- loadGenome()

server <- function(input, output, session) {
  filtered_phenotype <- filterSamplesServer("filter_samples", phenotype)
  filtered_expression <- filterTranscriptsServer("filter_transcripts", circ_cpm, filtered_phenotype)
  dimredServer("dimred", filtered_expression, filtered_phenotype)
  statisticsServer("statistics", filtered_expression, filtered_phenotype)
}

shinyApp(ui = ui, server = server)
