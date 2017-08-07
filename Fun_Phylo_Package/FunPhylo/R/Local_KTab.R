
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
