library(shiny)

# Set port
options(shiny.host = "0.0.0.0")
options(shiny.port = 8080)
runApp("src")
