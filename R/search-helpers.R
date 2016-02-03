# Functions for searching using the AISearch package

#' Create the initial search node
initSearchNode = function ( lll, reporters = getReporters( lll ) ){
  n = howManyActors( lll ) # Number of actors
  nR = length( reporters ) # Number of reporters
  
  the.dim.names = list( getActors( lll ), c( getActors( lll ), reporters ) )
  
  uncertain = ( 1 != diag( nrow = n, ncol = n + nR ) )
  dimnames( uncertain ) = the.dim.names
  adjacency = matrix( FALSE, nrow = n, ncol = n + nR, dimnames = the.dim.names )
  
  makeSearchNode ( lll, adjacency, ancestry )
}

makeSearchNode = function ( lll, adjacency, ancestry ){
  # Build uncertainty-based ancestry matrices
  possible.ancestors = adjacencyToAncestry( adjacency | uncertain )
  possible.nonancestors = !adjacencyToAncestry( adjacency & !uncertain )
  
  list(
    adjacency = adjacency,
    uncertain = uncertain,
    bound.vector = getScoreBounds( possible.ancestors, possible.nonancestors )
  )
}

#' We reach a search goal when nothing is uncertain, or the upperbound is the lowerbound
isGoalNode = function ( search.node ) {
  ( !any( search.node$uncertain ) ) || all( search.node$bound.vector["upper",] == search.node$bound.vector["lower",] )
}

#' Extract the cost bounds from the node
#' 
#' We negate the scores because the search algorithm uses costs while we score based on likelihoods
#' @param search.node The search node object
#' @return Cost bounds 
getCost = function ( search.node ) {
  c( "upper" = -sum( search.node$bound.vector["lower"] ),
     "lower" = -sum( search.node$bound.vector["upper"] ) )
}

#' Get children in the search graph
#' 
#' @param lll the LocalLogLikelihoods object
#' @param search.node the search node whose children we're generating
#' @return a list of search nodes (the children)
getChildNodes = function ( lll, search.node ) {
  
  n = howManyActors( lll )
  
  result = list()
  
  uncertain = search.node$uncertain
  ancestry = search.node$ancestry
  
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
    ancestry[ edge.row, edge.col ] = edge
    result = c( result, list( makeSearchNode( lll=lll, ancestry=ancestry, uncertain=uncertain ) ) )
  }
  
  result
}