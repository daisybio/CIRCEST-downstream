library(shiny)

aboutUI <- tabPanel(
  "About",
  titlePanel("About circRNA investigator"),
  mainPanel(
    h2("About circRNA investigator"),
    p("This is a Shiny application for investigating circRNA.")
  )
)
