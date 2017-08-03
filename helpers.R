# helper functions for Shiny_Phylo_Fun!

mod_func_phy_dist<-function (FDist, PDist, phyloWeight, p, ...) 
{
  if (phyloWeight < 0 | phyloWeight > 1) 
    stop("'phyloWeight' must be between 0 and 1")
  if (!is.numeric(p)) 
    stop("'p' must be a numeric")
  
  FDist <- FDist/max(FDist)
  PDist <- PDist/max(PDist)
  (phyloWeight * PDist^p + (1 - phyloWeight) * FDist^p)^(1/p)
}

make_local_phylo_dist <- function(focal.species, community.data,
                             phylo, traits){
  local.com <- dplyr::filter(community.data, 
                             exotic_species == focal.species) %>%
               .$community %>% as.character() %>%
               .[. %in% traits$Species.Name]
  
  out <- ape::drop.tip(phylo, setdiff(phylo$tip.label, local.com)) %>%
          cophenetic() %>% data.frame()
  
  return(out)
  
}

make_local_traits_ktab <- function(focal.species,community.data,
                                   trait.data, traits){
  
  wood <- trait.data[trait.data$Species.Name == focal.species, 'Woody']
  whole.com <- dplyr::filter(community.data,
                             exotic_species == focal.species) %>%
               .$community
  
  local.com <- trait.data[trait.data$Woody == wood &
                          trait.data$Species.Name %in% whole.com, ]
  
  if(wood == 0){
    ContNames <- HerbTraitNames %in% traits
  } else {
    ContNames <- WoodyTraitNames %in% traits
  }
  LHNames <- LHTNames %in% traits
  GFNames <- GFTNames %in% traits
  
  ### INCOMPLETE ###
}


