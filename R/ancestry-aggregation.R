## Functions having to do with ancestry aggregation etc.
##

#' Aggregate individual reporter-ancestries into one master ancestry
#' 
#' @param ancestries - a actors x actors+1 x reporters array of individual
#' reporter-ancestries.
#' @return a single ancestry (with possibly duplicated actors) which is a matrix of
#' size (actors) x (new actors + reporters)
aggregateAncestries = function ( ancestries ){
  
  # Get dimensions and names
  ds = dim( ancestries )
  nA = ds[1]
  nR = ds[3]
  
  dn = dimnames( ancestries )
  if ( ! is.null( dn ) ){
    actors = dn[[ 1 ]]
    reporters = dn[[ 3 ]]
    actorMap = vector( mode = "character", length = 0 )
  } else {
    actors = 1:nA
    reporters = 1:nR
    actorMap = vector( mode = "numeric", length = 0 )
  }  
  
  actor.ancestry = reporter.ancestry = matrix( nrow = nA, ncol = 0 )
  rownames( actor.ancestry ) = colnames( reporter.ancestry ) = actors
  
  for ( reporter in reporters ) {
    ancestry = ancestries[, nA+1, reporter ]
    for ( actor in actors( ancestry ) ) {
      #
    }
  }
}