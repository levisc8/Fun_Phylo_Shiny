
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
