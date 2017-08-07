make_local_trait_dist <- function(focal.species, community.data, trait.data, traits, scale){
  

  ktab <- make_local_traits_ktab(focal.species, community.data, trait.data, 
                                 traits)
  
  out <- ade4::dist.ktab(ktab$KTab, type = ktab$VarTypes,
                         option = scale)
  
  return(out)
  
  
}
