library(shiny)

correlationUI <- card(
  card_header("Correlation options"),
  card_body(
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
  )
)

differentialUI <- card(
  card_header("Differential expression options"),
  card_body(
    "Hello world",
    input_task_button("run_differential", "Run differential")
  )
)

statisticsUI <- tabPanel(
  "Statistics",
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "test_type",
        "Test type:",
        choices = c("correlation", "differential"),
        selected = "correlation"
      ),
      uiOutput("statisticsUI")
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
          ),
          downloadButton("download_correlation", "Download results")
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
