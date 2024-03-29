library(shiny)
library(Hmisc) # for rcorr
library(MASS) # for mvrnorm
library(ggplot2) # for plotting
library(ggpubr) # for ggarrange
library(fGarch) # for skewed distribution
library(WRS2) # for wincor and pbcor
library(correlation) # for multilevel correlations
source("SpurCorr.R")

# User interface ----
ui <- fluidPage(
  titlePanel("Simulate Correlation"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Create your own (spurious) correlation(s)"),
      
      selectInput("var", 
                  label = "Choose a data type to display",
                  choices = c("Correlation", "Correlation with outlier",
                              "Correlation with subgroups"),
                  selected = "Correlation"),
      selectInput("Solution", 
                  label = "Account for data irregularity",
                  choices = c("no", "outlier: robust correlation","subgroups: regression","subgroups: multilevel correlation"),
                  selected = "no"),
      numericInput("Nsim", 
                  label = "Number of simulations",
                  min = 1, max = 10000, value = 100),
      numericInput("Rvalue", 
                label = "Actual Correlation",
                min = -1, max = 1, value = 0),
      numericInput("PopSize", 
                label = "Population size",
                min = 2, max = 1000, value = 15),
    sliderInput("SDdist", 
                label = "distance from mean: outlier/Population",
                min = 0, max = 5, value = 3),
    actionButton("action", label = "Generate again"),
    selectInput("BinSize", "Bin Size (for p-values):",
                c("0.05" = "wide",
                  "0.01" = "narrow"))
  ),  
  mainPanel(plotOutput("map"))
  )
)

# Server logic ----
server <- function(input, output) {
  output$map <- renderPlot({
    input$action
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
                        Solution=input$Solution,
                        BS=input$BinSize))
  })
}

# Run app ----
shinyApp(ui, server)