
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("Invasive Plant Species and the Role of Competition and Novelty"),


  # Sidebar with a slider input for number of bins
  fluidRow(
      column(3,
             radioButtons("scale", label = "Spatial Scale",
                          choices = list("Local (Plot Level)" = 'loc',
                                         "Regional (County Level)" = 'reg')),
             
             radioButtons("met.phylo", label = "Phylogenetic Metric",
                          choiceNames = c("Nearest Neighbor Distance",
                                         "Mean Pairwise Distance"),
                          choiceValues = c('NND', 'MPD')),
             
             radioButtons("met.inv", label = "Metric of Invasiveness",
                          choices = list("Per-capita Growth Rates" = 'lambda',
                                         "Expert Classification" = 'mepp')),
             
             sliderInput("Little.a", "Phylogenetic Scaling Parameter",
                         min = 0, max = 1, value = .5, step = .05)),
             
             helpText('This app employs the methods described',
                      ' in Cadotte et al 2013 and Thorn et al 2015',
                      ' to demonstrate how the relationship between',
                      ' novelty and invasiveness is mediated by competition.',
                      ' The goal here is to create an interactive',
                      ' version of that analysis that allows a user ',
                      'to explore how relationships can change when ',
                      'incorporating information ',
                      'on both phylogeny and functional traits. Future directions',
                      'include generalizing this app so other users ',
                      'can put in their own data and phylogeny and',
                      'conduct analyses using this interactive platform.')
    ),
    
    column(3,
        checkboxGroupInput("traits", label = "Functional Traits",
                           choices = list('SLA' = "SLA", 'Height' = "Height",
                                          'Tough' = "Tough",
                                          'Wood Density (if applicable)' = "WoodDens", 
                                          'Flower Period' = "Flower.Period", 
                                          'Growth Form' = 'Growth_Form',
                                          'Dispersal Mechanism' = 'Disp_Mech'),
                           selected = c('SLA', 'Height', 'Tough','Flower.Period'))

    ),
    mainPanel(
      plotOutput('figure1')#,
      # tableOutput('traits')
    )
    
  )
)
