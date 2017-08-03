
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("Functional-Phylogenetic Distances"),

  # Sidebar with a slider input for number of bins
  fluidRow(
      column(3,
             radioButtons("scale", label = "Spatial Scale",
                          choices = list("Local (Plot Level)",
                                         "Regional (County Level)")),
             
             radioButtons("met.phylo", label = "Phylogenetic Metric",
                          choices = list("Nearest Neighbor Distance",
                                         "Mean Pairwise Distance")),
             
             radioButtons("met.inv", label = "Metric of Invasiveness",
                          choices = list("Per-capita Growth Rates",
                                         "Expert Classification")),
             
             sliderInput("Little.a", "Phylogenetic Scaling Parameter",
                         min = 0, max = 1, value = 0.5, step = .05)),
             
             helpText('This analysis follows the methods described',
                      ' in Cadotte et al 2013 and Thorn et al 2015.',
                      ' The goal here is to create an interactive',
                      ' version of that analysis that allows a user ',
                      'to explore how previously unexposed relationships',
                      ' can become important when incorporating information ',
                      'on both phylogeny and functional traits.')
    ),
    
    column(3,
        checkboxGroupInput("traits", label = "Functional Traits",
                           choices = list('SLA', 'Height', 'Tough',
                                          'Wood Density', 'Woody',
                                          'Flower Period', 'Clonal',
                                          'Legume', 'Growth Form',
                                          'Dispersal Mechanism'),
                           selected = c('SLA', 'Height', 'Tough',
                                        'Wood Density', 'Woody',
                                        'Flower Period', 'Clonal',
                                        'Legume', 'Growth Form',
                                        'Dispersal Mechanism'))

    ),
    mainPanel(
      plotOutput('figure1')
    )
    
  )
)
