library(shiny)

aboutUI <- tabPanel(
  "About",
  titlePanel("About circRNA investigator"),
  mainPanel(
    p(paste(
      "This is a Shiny application for investigating",
      "the output of the nf-core/circrna pipeline."
    )),
  )
)
