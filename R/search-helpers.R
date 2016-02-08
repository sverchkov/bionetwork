# Functions for searching using the AISearch package

#' Create the initial search node
#' 
#' Creates the initial (all-uncertain) node.
#' @param lll LocalLogLikelihoods object
#' @param reporters array of reporter names to use (optional, default pulled from lll)
#' @return a search node (specified by adjacency and uncertainty matrices and a vector
#' of reporter contributions to the score)
#' @export
initSearchNode = function ( lll, reporters = getReporters( lll ) ){
  n = howManyActors( lll ) # Number of actors
  nR = length( reporters ) # Number of reporters
  
  the.dim.names = list( getActors( lll ), c( getActors( lll ), reporters ) )
  
  uncertain = ( 1 != diag( nrow = n, ncol = n + nR ) )
  dimnames( uncertain ) = the.dim.names
  adjacency = matrix( FALSE, nrow = n, ncol = n + nR, dimnames = the.dim.names )
  
  makeSearchNode ( lll, adjacency, uncertain )
}

#' Makes a search node
#' 
#' @param lll LocalLogLikelihoods object
#' @param adjacency adjacency matrix
#' @param uncertain uncertainty matrix
#' @return a search node (list made up of adjacency, uncertainty, and per-reporter
#' score vector)
makeSearchNode = function ( lll, adjacency, uncertain ){
  # Build uncertainty-based ancestry matrices
  possible.ancestors = adjacencyToAncestry( adjacency | uncertain )
  possible.nonancestors = !adjacencyToAncestry( adjacency & !uncertain )
  
  list(
    adjacency = adjacency,
    uncertain = uncertain,
    bound.vector = getScoreBounds( lll, possible.ancestors, possible.nonancestors )
  )
}

#' Goal node check
#' 
#' We reach a search goal when nothing is uncertain, or the upperbound is the lowerbound
#' @param search.node the search node
#' @return TRUE iff this is a goal node
#' @export
isGoalNode = function ( search.node ) {
  ( !any( search.node$uncertain ) ) || all( search.node$bound.vector["upper",] == search.node$bound.vector["lower",] )
}

#' Extract the cost bounds from the node
#' 
#' We negate the scores because the search algorithm uses costs while we score based on likelihoods
#' @param search.node The search node object
#' @return Cost bounds
#' @export
getCost = function ( search.node ) {
  c( "upper" = -sum( search.node$bound.vector["lower", ] ),
     "lower" = -sum( search.node$bound.vector["upper", ] ) )
}

#' Get children in the search graph
#' 
#' @param lll the LocalLogLikelihoods object
#' @param search.node the search node whose children we're generating
#' @return a list of search nodes (the children)
#' @export
getChildNodes = function ( lll, search.node ) {
  
  n = howManyActors( lll )
  
  result = list()
  
  uncertain = search.node$uncertain
  adjacency = search.node$adjacency
  
  # Figure out which node to assign
  edge.col =
    if ( any( candidates <- apply( uncertain[ , -( 1:n ) ], 2, any ) ) ){
      variances = search.node$bound.vector["upper", candidates] - search.node$bound.vector["lower", candidates]
      names( which.max( variances ) )
    } else {
      names( which.max( apply( uncertain, 2, any ) ) )
    }
  edge.row = which.max( uncertain[ , edge.col ] )
  
  uncertain[ edge.row, edge.col ] = FALSE
  
  for ( edge in c( FALSE, TRUE ) ){
    adjacency[ edge.row, edge.col ] = edge
    result = c( result, list( makeSearchNode( lll=lll, adjacency=adjacency, uncertain=uncertain ) ) )
  }
  
  result
}