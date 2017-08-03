# tests for server

rm(list = ls(all = T))
library(shiny)
library(tidyverse)
library(magrittr)
library(ade4)
library(ape)
library(pez)
source('helpers.R')

load('data/data_list.RData')

communities <- data$communities
phylo <- data$phylo
tyson <- data$spp.list
demo <- data$demo.data
trait.data <- data$traits

traits <- names(trait.data[-1])

for(i in 1:10){
  trait.test <- traits[base::sample(1:24,9)]
  
  for(x in unique(demo$Species)){
    cat(x,"    ", i,'\n')
    test.local <- make_local_traits_ktab(x, community.data = communities,
                                         trait.data = trait.data, traits = traits)
  }
  
  test.regional <- make_regional_traits_ktab(trait.data, traits)
  
  
}
