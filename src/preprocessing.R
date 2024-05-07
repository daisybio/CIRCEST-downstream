library(shiny)

preprocessingUI <- tabPanel(
  "Preprocessing",
  sidebarLayout(
    sidebarPanel(
      card(
        card_header("Filter transcripts"),
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
          ),
          radioButtons("transcript_types", "Transcript types",
            choices = c("circular", "linear", "both"),
            selected = "circular"
          ),
          textOutput("filtered_description")
        )
      )
    ),
    mainPanel(
      card(
        card_header("Settings"),
        card_body(
          uiOutput("colorings")
        )
      ),
      card(
        card_header("PCA"),
        card_body(
          plotlyOutput("plotPCA"),
        )
      ),
      card(
        card_header("UMAP"),
        card_body(
          plotlyOutput("plotUMAP")
        )
      )
    )
  )
)
