library(shiny)

correlationUI <- tabPanel(
  "Correlation",
  sidebarLayout(
    sidebarPanel(
      selectizeInput("cor_x",
        "Correlation X:",
        choices = NULL
      ),
      selectInput("cor_type",
        "Correlation type:",
        choices = c("pearson", "spearman"),
        selected = "pearson"
      ),
      input_task_button("run_correlation", "Run correlation")
    ),
    mainPanel(
      sliderInput("cor_alpha",
        "Alpha",
        min = 0,
        max = 1,
        value = 0.05
      ),
      plotlyOutput("cor_volcano"),
      plotlyOutput("cor_heatmap")
    )
  )
)
