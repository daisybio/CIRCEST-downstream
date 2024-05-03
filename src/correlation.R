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
      verbatimTextOutput("summary"),
      plotlyOutput("cor_volcano")
    )
  )
)
