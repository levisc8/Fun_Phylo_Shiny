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
  
  WoodyTraitNames <- c('SLA','Tough','WoodDens')
  HerbTraitNames <- c("Height","SLA","Tough")

  GrowthForm <- c("Stemmed_Herb",
                  "Tree", "Rosette",
                  "Vine", "SubShrub",
                  "Shrub",
                  "Elongated_Leafy_Rhizomatous",
                  "N_Fixer")
  
  DispersalTraitNames <- c("Subterranean",
                           "Unassisted",
                           "Wind",
                           "ExoZoochory",
                           "EndoZoochory",
                           "Ballistic",
                           "Hoarding",
                           "Myrmecochory",               
                           "Water",
                           "Clonal")
  
  wood <- trait.data[trait.data$Species.Name == focal.species, 'Woody']
  whole.com <- dplyr::filter(community.data,
                             exotic_species == focal.species) %>%
               .$community
  
  local.com <- trait.data[trait.data$Woody == wood &
                          trait.data$Species.Name %in% whole.com, ]
  
  if(wood == 0){
    ContNames <- HerbTraitNames[HerbTraitNames %in% traits]
  } else {
    ContNames <- WoodyTraitNames[WoodyTraitNames %in% traits]
  }
  GFNames <- GrowthForm[GrowthForm %in% traits]
  DispNames <- DispersalTraitNames[DispersalTraitNames %in% traits]
  
  if(length(ContNames) > 0){
    ContTraits <- data.frame(local.com[ ,ContNames])
  }
  if(length(GFNames) > 0){
    GFTraits <- data.frame(local.com[ ,GFNames])
  } 
  if(length(DispNames) > 0){
    DispTraits <- data.frame(local.com[ ,DispNames])
  }
  if("Flower.Period" %in% traits){
    Flower.Period <- data.frame(local.com[ ,'Flower.Period'])
  }
  
  VarTypes <- NULL
  
  
  if(exists("ContTraits")){
    VarTypes <- c(VarTypes, "Q")
  } else {
    ContTraits <- NULL
  }
  
  if(exists("GFTraits")){
    GFKTab <- ade4::prep.binary(GFTraits,
                                col.blocks = dim(GFTraits)[2],
                                label = "Growth_Form")
    VarTypes <- c(VarTypes, "B")
  } else {
    GFKTab <- NULL
  }
  
  if(exists("DispTraits")){
    DispKTab <- ade4::prep.binary(DispTraits, 
                                  col.blocks = dim(DispTraits)[2],
                                  label="Dispersal_Traits")
    VarTypes <- c(VarTypes, "B")
    
  } else {
    DispKTab <- NULL
  }
  
  if(exists("Flower.Period")){
    FlowerKTab <- ade4::prep.circular(Flower.Period,
                                      rangemin=1, rangemax=12)
    VarTypes <- c(VarTypes, "C")
    
  } else {
    FlowerKTab <- NULL
  }

  
  forKTab <- check_empty_elements(ContTraits, GFKTab,
                                   DispKTab, FlowerKTab)

  KTab <- ade4::ktab.list.df(forKTab,
                             rownames = local.com$Species.Name)
  
  out <- list(KTab = KTab, VarTypes = VarTypes)
  
  return(out)
  
}

check_empty_elements <- function(...){
  rmlist <- NULL
  iniList <- list(...)
  
  for(i in 1:length(iniList)){
    if(is.null(iniList[[i]])){
      rmlist <- c(rmlist, i)
    }
  }
  
  if(!is.null(rmlist)){
    outList <- iniList[-c(rmlist)]
  } else {
    outList <- iniList
  }
  
  return(outList)
}

make_regional_traits_ktab <- function(traits.data, traits){
  
  ContTraitNames <- c('SLA','Tough','WoodDens','Height')
  
  GrowthForm <- c("Woody", "Stemmed_Herb",
                  "Tree", "Rosette",
                  "Vine", "SubShrub",
                  "Shrub",
                  "Elongated_Leafy_Rhizomatous",
                  "N_Fixer")
  
  DispersalTraitNames <- c("Subterranean",
                           "Unassisted",
                           "Wind",
                           "ExoZoochory",
                           "EndoZoochory",
                           "Ballistic",
                           "Hoarding",
                           "Myrmecochory",               
                           "Water",
                           "Clonal")

  spp.names <- trait.data$Species.Name
  
  ContNames <- ContTraitNames[ContTraitNames %in% traits]
  GFNames <- GrowthForm[GrowthForm %in% traits]
  DispNames <- DispersalTraitNames[DispersalTraitNames %in% traits]
  
  if(length(ContNames) > 0){
    ContTraits <- data.frame(trait.data[ ,ContNames])
  }
  if(length(GFNames) > 0){
    GFTraits <- data.frame(trait.data[ ,GFNames])
  } 
  if(length(DispNames) > 0){
    DispTraits <- data.frame(trait.data[ ,DispNames])
  }
  if("Flower.Period" %in% traits){
    Flower.Period <- data.frame(trait.data[ ,'Flower.Period'])
  }
  
  VarTypes <- NULL
  
  if(exists("ContTraits")){
    VarTypes <- c(VarTypes, "Q")
  } else {
    ContTraits <- NULL
  }
  
  if(exists("GFTraits")){
    GFKTab <- ade4::prep.binary(GFTraits,
                                col.blocks = dim(GFTraits)[2],
                                label = "Growth_Form")
    VarTypes <- c(VarTypes, "B")
  } else {
    GFKTab <- NULL
  }
  
  if(exists("DispTraits")){
    DispKTab <- ade4::prep.binary(DispTraits, 
                                  col.blocks = dim(DispTraits)[2],
                                  label="Dispersal_Traits")
    VarTypes <- c(VarTypes, "B")
    
  } else {
    DispKTab <- NULL
  }
  
  if(exists("Flower.Period")){
    FlowerKTab <- ade4::prep.circular(Flower.Period,
                                      rangemin=1, rangemax=12)
    VarTypes <- c(VarTypes, "C")
    
  } else {
    FlowerKTab <- NULL
  }
  
  
  forKTab <- check_empty_elements(ContTraits, GFKTab,
                                  DispKTab, FlowerKTab)
  
  KTab <- ade4::ktab.list.df(forKTab,
                             rownames = spp.names)
  
  out <- list(KTab = KTab, VarTypes = VarTypes)
  
  return(out)
  
}

rarefy_FPD <- function(focal.species, phylo.mat, fun.mat,
                               metric = c("MPD", "NND"), 
                               n.resamp = 1000, n.rare) {
  mpd.tf <- "MPD" %in% metric
  nnd.tf <- "NND" %in% metric
  
  if(mpd.tf){
    rare.mpd <- rep(NA, n.resamp)
  }
  if(nnd.tf){
    rare.bl <- rep(NA, n.resamp)
  }
  
  diag(dist.mat) <- NA
  focal.column <- dist.mat[ ,focal.species]
  
  for(i in 1:n.resamp){
    resamp.x <- base::sample(1:length(focal.column),
                             size = n.rare,
                             replace = FALSE)
    
    if(mpd.tf) {
      rare.mpd[i] <- mean(focal.column[resamp.x], na.rm = TRUE)
    }
    if(nnd.tf){
      rare.bl[i] <- min(focal.column[resamp.x], na.rm = TRUE)
    }
  }
  
  out <- list()
  if(mpd.tf){
    out$rare.mpd <- mean(rare.mpd)
  }
  if(nnd.tf){
    out$rare.nnd <- mean(rare.bl)
  }
  
  return(out)
  
}











