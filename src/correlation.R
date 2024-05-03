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
      card(
        card_header("Plot options"),
        card_body(
          sliderInput("cor_alpha",
            "Alpha",
            min = 0,
            max = 1,
            value = 0.05
          )
        )
      ),
      card(
        card_header("Volcano plot"),
        card_body(
          plotlyOutput("cor_volcano")
        )
      ),
      card(
        card_header("Heatmap"),
        card_body(
          plotlyOutput("cor_heatmap")
        )
      )
    )
  )
)
