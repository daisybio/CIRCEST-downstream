library(shiny)

datasets <- list.dirs(path = "../data", full.names = TRUE, recursive = FALSE)

preprocessingUI <- tabPanel(
  "Preprocessing",
  titlePanel("circRNA investigator"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("bins",
        "Number of bins:",
        min = 1,
        max = 50,
        value = 30
      ),
      selectInput("dataset",
        "Select dataset:",
        choices = datasets
      )
    ),
    mainPanel(
      plotOutput("distPlot"),
      textOutput("datadescription")
    )
  )
)
