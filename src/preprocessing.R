library(shiny)

datasets <- list.dirs(path = "../data", full.names = TRUE, recursive = FALSE)

preprocessingUI <- tabPanel(
  "Preprocessing",
  sidebarLayout(
    sidebarPanel(
      card(
        card_header("Select dataset"),
        card_body(
          selectInput("dataset",
            "Select dataset:",
            choices = datasets,
            selected = "../data/test"
          )
        )
      ),
      card(
        card_header("Filter circRNAs"),
        card_body(
          sliderInput("min_count",
            "Min count",
            min = 0,
            max = 100,
            value = 50
          ),
          sliderInput("min_samples_pct",
            "Min samples %",
            min = 0,
            max = 100,
            value = 20.0
          )
        )
      )
    ),
    mainPanel(
      uiOutput("colorings"),
      plotlyOutput("plotPCA"),
      plotlyOutput("plotUMAP")
    )
  )
)
