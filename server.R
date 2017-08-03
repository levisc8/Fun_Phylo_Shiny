
# code that has to run regardless of user goes here
# e.g. data assignment (though that may change if
# generalized for other users). 

library(shiny)
library(tidyverse)
library(magrittr)
library(ade4)
library(ape)
library(pez)

load('data/data_list.RData')

communities <- data$communities
phylo <- data$phylo
tyson <- data$spp.list
demo <- data$demo.data
traits <- data$traits

shinyServer(function(input, output) {
  
  # some code for plot goes here
  
  renderPlot({
    
    # most code for plot goes here
    
    output$figure1 <- FigureForOutput})
})
