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
# test trait ktab functions
for(i in 1:10){
  trait.test <- traits[base::sample(1:24,9)]
  
  for(x in unique(demo$Species)){
    cat(x,"    ", i,'\n')
    test.local <- make_local_traits_ktab(x, community.data = communities,
                                         trait.data = trait.data, traits = traits)
  }
  
  test.regional <- make_regional_traits_ktab(trait.data, traits)
  
  
}


# test rarefying functions

local.com <- filter(communities, exotic_species == 'Ailanthus_altissima')

local.phy <- drop.tip(phylo, setdiff(phylo$tip.label, local.com$community))
test.dist <- data.frame(cophenetic(local.phy))
diag(test.dist) <- NA
min.thresh <- min(test.dist$Ailanthus_altissima, na.rm = T)

n.runs <- 200
for(i in 1:n.runs){
  test <- rarefy_phylo_dists('Ailanthus_altissima', test.dist, n.rare = 15)
  
  if(test$rare.nnd < min.thresh) warning('You fucked up')
}