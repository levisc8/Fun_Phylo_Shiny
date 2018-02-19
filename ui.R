
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
    column(2,
           selectInput("plot", "Plot Type", 
                       choices = list('GLM' = 'FS',
                                      'Explanatory Power' = 'lil.a'),
                       selectize = FALSE)),
    column(2,     
           radioButtons("scale", label = "Spatial Scale",
                        choices = list("Local (Plot Level)" = 'loc',
                                       "Regional (County Level)" = 'reg'),
                        selected = 'loc')),
    column(4,
           checkboxGroupInput("traits", label = "Functional Traits",
                              choices = list('SLA' = "SLA", 'Height' = "Height",
                                             'Leaf Toughness' = "Tough",
                                             'Wood Density (if applicable)' = "WoodDens", 
                                             'First Flowering' = "Flower.Period", 
                                             'Growth Form' = 'Growth_Form',
                                             'Dispersal Mechanism' = 'Disp_Mech'),
                              selected = c('SLA', 'Height', 
                                           'Tough','Flower.Period'),
                              inline = TRUE)),
    column(2, 
           checkboxInput('AW', label = 'Weight Local Species Abundances',
                         value = FALSE)),
    column(2,
           checkboxInput('log', label = 'Log transform abundance weighted data?'),
           value = FALSE)
    
    
  ),

  # Sidebar with a slider input for number of bins
  conditionalPanel(condition = "input.plot == 'FS'",
    fluidRow(
      column(3,
             radioButtons("met.phylo", label = "Phylogenetic Metric",
                          choiceNames = c("Nearest Neighbor Distance",
                                          "Mean Pairwise Distance"),
                          choiceValues = c('NND', 'MPD'))),
      column(3,
             radioButtons("met.inv", label = "Metric of Invasiveness",
                          choices = list("Per-capita Growth Rates" = 'lambda',
                                         "Expert Classification" = 'mepp'))),
      column(3,      
             sliderInput("Little.a", "Phylogenetic Scaling Parameter",
                         min = 0, max = 1, value = .5, step = .05))
      
    )
    
  ),
  conditionalPanel(condition = "input.plot == 'lil.a'",
   fluidRow(
     column(3,
            sliderInput("res", label = 'Breaks for Phylogenetic Scaling Parameter',
                        min = 5, max = 41, value = 11, step = 2)

            
     )
     
   )
    
  ),
    mainPanel(
      plotOutput('figure1')#,
      # tableOutput('table1')
    )
    
  )
)
