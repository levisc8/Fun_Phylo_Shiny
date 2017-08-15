
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(markdown)

shinyUI(fluidPage(

  # Application title
  titlePanel("Invasive Plant Species and the Role of Competition and Novelty"),
  includeMarkdown('Readme_general.md'),
  fluidRow(
    column(3,
           selectInput("plot", "Plot Type", 
                       choices = list('Focal Species' = 'FS',
                                      'R^2 ~ a' = 'lil.a'),
                       selectize = FALSE),
           
           radioButtons("scale", label = "Spatial Scale",
                        choices = list("Local (Plot Level)" = 'loc',
                                       "Regional (County Level)" = 'reg'),
                        selected = 'loc'),
           checkboxGroupInput("traits", label = "Functional Traits",
                              choices = list('SLA' = "SLA", 'Height' = "Height",
                                             'Tough' = "Tough",
                                             'Wood Density (if applicable)' = "WoodDens", 
                                             'Flower Period' = "Flower.Period", 
                                             'Growth Form' = 'Growth_Form',
                                             'Dispersal Mechanism' = 'Disp_Mech'),
                              selected = c('SLA', 'Height', 
                                           'Tough','Flower.Period')))
  ),
  # Sidebar with a slider input for number of bins
  conditionalPanel(condition = "input.plot == 'FS'",
    fluidRow(
      column(3,
             radioButtons("met.phylo", label = "Phylogenetic Metric",
                          choiceNames = c("Nearest Neighbor Distance",
                                          "Mean Pairwise Distance"),
                          choiceValues = c('NND', 'MPD')),
             
             radioButtons("met.inv", label = "Metric of Invasiveness",
                          choices = list("Per-capita Growth Rates" = 'lambda',
                                         "Expert Classification" = 'mepp')),
             
             sliderInput("Little.a", "Phylogenetic Scaling Parameter",
                         min = 0, max = 1, value = .5, step = .05))
      
    )
    
  ),
  conditionalPanel(condition = "input.plot == 'lil.a'",
   fluidRow(
     column(3,
            sliderInput("res", label = 'Breaks for Phylogenetic Scaling Parameter',
                        min = 5, max = 41, value = 21, step = 2)

            
     )
     
   )
    
  ),
    mainPanel(
      plotOutput('figure1')#,
      # textOutput('traits')
    )
    
  )
)
