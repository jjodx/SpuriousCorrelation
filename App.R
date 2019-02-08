library(shiny)
library(Hmisc) # for rcorr
library(MASS) # for mvrnorm
library(ggplot2) # for plotting
library(ggpubr) # for ggarrange
library(fGarch) # for skewed distribution
library(WRS2) # for wincor and pbcor
source("SpurCorr.R")

# User interface ----
ui <- fluidPage(
  titlePanel("Simulate Correlation"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Create your own (spurious) correlation(s)"),
      
      selectInput("var", 
                  label = "Choose a variable to display",
                  choices = c("Correlation", "Correlation with outlier",
                              "Correlation with subgroups"),
                  selected = "Correlation"),
      selectInput("Solution", 
                  label = "Use robust correlations?",
                  choices = c("yes", "no"),
                  selected = "no"),
      sliderInput("Nsim", 
                  label = "Number of simulations",
                  min = 1, max = 10000, value = 500),
    sliderInput("Rvalue", 
                label = "Actual Correlation",
                min = 0, max = 1, value = 0),
    sliderInput("PopSize", 
                label = "Number of participants",
                min = 2, max = 100, value = 20),
    sliderInput("SDdist", 
                label = "distance from mean: outlier/Population",
                min = 1, max = 5, value = 3)
  ),  
  mainPanel(plotOutput("map"))
  )
)

# Server logic ----
server <- function(input, output) {
  output$map <- renderPlot({
    PR <- switch(input$var, 
                   "Correlation" = "none",
                   "Correlation with outlier"= "outlier",
                   "Correlation with subgroups" = "subpop"
                   )
 
    report(SimulateCorr(PopSize=input$PopSize,
                        Nsim=input$Nsim,
                        Rval=input$Rvalue,
                        Problem=PR,
                        SDdist=input$SDdist,
                        Skewness=2,
                        Solution=input$Solution))
  })
}

# Run app ----
shinyApp(ui, server)