library(shiny)
library(bslib)
library(DESeq2)

source("preprocessing.R")
source("about.R")

# Define UI for app that draws a histogram ----
ui <- navbarPage(
  "circRNA investigator",
  preprocessingUI,
  aboutUI
)

normalizeSE <- function(se) {
  colData <- colData(se)

  dds <- DESeqDataSet(se, )
}

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  output$distPlot <- renderPlot({
    x <- faithful$waiting
    bins <- seq(min(x), max(x), length.out = input$bins + 1)

    hist(x,
      breaks = bins, col = "#007bc2", border = "white",
      xlab = "Waiting time to next eruption (in mins)",
      main = "Histogram of waiting times"
    )
  })

  dataset <- reactive({
    se <- readRDS(paste0(input$dataset, "/tx.rds"))
    rownames(colData(se)) <- colData(se)$names
    colData(se)$names <- NULL
    se
  })

  normalized <- reactive({
    se <- dataset()

    # Keep only columns with more than on unique value
    colDataFiltered <- colData(se)[, sapply(colData(se), function(x) length(unique(x)) > 1)]
    names <- colnames(colDataFiltered)
    design <- formula(paste0("~", paste(names, collapse = "+")))

    # Normalize
    dds <- DESeqDataSetFromMatrix(
      round(assay(se), 0),
      colDataFiltered,
      design = design
    )
    dds <- DESeq(dds)

    assay(se, "norm") <- counts(dds, normalized = TRUE)

    se
  })

  output$datadescription <- renderPrint({
    normalized()
  })
}

shinyApp(ui = ui, server = server)
