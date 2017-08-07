#' @title Rarefied Function-Phylogenetic Distances
#' @description 
#' \code{rare_FPD} resamples a phylogenetic and functional distance matrix
#' a set number of times so that communities of different sizes can be 
#' compared to each other.
#' 
#' @param focal.species The name of the focal species in the 
#' community to be rarefied.
#' @param phylo.mat,fun.mat A phylogenetic or functional distance
#' matrix with species names as row names and column names. 
#' Can either be a \code{data frame}, \code{matrix}, or \code{dist}.
#' @param metric Either MPD (\code{MPD}) or nearest neighbor distance
#' (\code{NND}). The default is to calculate both metrics.
#' @param n.resamp The number of rarefied resamplings to perform. 
#' Default is 1000. 
#' @param n.rare The number of species to rarefy the community down to.
#' Must be less than the total in the community.
#' @param a The phylogenetic scaling factor for the calculation of
#' functional-phylogenetic distances
#'  
#' @export
#' 
rarefy_FPD <- function(focal.species, phylo.mat, fun.mat,
                       metric = c("MPD", "NND"),
                       n.resamp = 1000, n.rare, 
                       a) {
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

