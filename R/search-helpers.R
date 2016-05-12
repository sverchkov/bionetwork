# Functions for searching using the AISearch package

#' Create the initial search node
#' 
#' Creates the initial (all-uncertain) node.
#' @param lll LocalLogLikelihoods object
#' @param reporters array of reporter names to use (optional, default pulled from lll)
#' @return a search node (specified by adjacency and uncertainty matrices and a vector
#' of reporter contributions to the score)
#' @export
initSearchNode = function ( lll, reporters = getReporters( lll ), use.sparse = TRUE ){
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
#' @param auto.sparse whether to try and use sparse matrices
#' @return a search node (list made up of adjacency, uncertainty, and per-reporter
#' score vector)
#' @export
makeSearchNode = function ( lll, adjacency, uncertain, auto.sparse = TRUE ){
  # Build uncertainty-based ancestry matrices
  possible.ancestors = adjacencyToAncestry( adjacency | uncertain )
  possible.nonancestors = !adjacencyToAncestry( adjacency & !uncertain )
  
  bounds = getScoreBounds( lll, possible.ancestors, possible.nonancestors )
  
  if ( auto.sparse ){
    sparse.adj = ( sum( adjacency )*10 > length( adjacency ) )
    sparse.unc = ( sum( uncertain )*10 > length( uncertain ) )
    list(
      adjacency = Matrix::Matrix( adjacency, sparse = sparse.adj ),
      uncertain = Matrix::Matrix( uncertain, sparse = sparse.unc ),
      bound.vector = bounds
    )
  } else
    list(
      adjacency = adjacency,
      uncertain = uncertain,
      bound.vector = bounds
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

#' Get children in the search graph (wide version)
#' 
#' @param lll the LocalLogLikelihoods object
#' @param search.node the search node whose children we're generating
#' @return a list of search nodes (the children)
#' @export
getChildNodesWide = function ( lll, search.node ) {
  c(
    lapply( which( search.node$uncertain ), function( index ) {
      uncertain = search.node$uncertain
      uncertain[ index ] = FALSE
      adjacency = search.node$adjacency
      adjacency[ index ] = FALSE
      makeSearchNode( lll=lll, adjacency=adjacency, uncertain=uncertain )
    } ),
    lapply( which( search.node$uncertain ), function( index ) {
      uncertain = search.node$uncertain
      uncertain[ index ] = FALSE
      adjacency = search.node$adjacency
      adjacency[ index ] = TRUE
      makesearchNode( lll=lll, adjacency=adjacency, uncertain=uncertain )
    } )
  )
}

#' Get children in the search graph
#' 
#' @param lll the LocalLogLikelihoods object
#' @param search.node the search node whose children we're generating
#' @param prioritize how to pick neighbors (by upperbound (default), or difference between upper and lower bound (variance) in the bound vector of the current node.)
#' @return a list of search nodes (the children)
#' @export
getChildNodes = function ( lll, search.node, prioritize = "upperbound" ) {
  
  n = howManyActors( lll )
  
  result = list()
  
  uncertain = search.node$uncertain
  adjacency = search.node$adjacency
  
  # Figure out which node to assign
  edge.col =
    if ( any( candidates <- apply( uncertain[ , -( 1:n ) ], 2, any ) ) ){
      priorities =
        if ( prioritize == "variance" )
          search.node$bound.vector["upper", candidates] - search.node$bound.vector["lower", candidates]
        else # upperbound
          search.node$bound.vector["upper", candidates]
      names( which.max( priorities ) )
    } else {
      names( which.max( apply( uncertain, 2, any ) ) )
    }
  
  if ( length( edge.col ) >= 1 ) {
    
    edge.row = which.max( uncertain[ , edge.col ] )
    
    uncertain[ edge.row, edge.col ] = FALSE
    
    for ( edge in c( FALSE, TRUE ) ){
      adjacency[ edge.row, edge.col ] = edge
      result = c( result, list( makeSearchNode( lll=lll, adjacency=adjacency, uncertain=uncertain ) ) )
    }
  }
  
  result
}

#' Simple node scorer
#' 
#' @param lll the LocalLogLikelihoods object
#' @param adjacency the adjacency matrix
#' @export
scoreSimply = function( lll, adjacency ){
  # Build ancestry matrices
  ancestors = adjacencyToAncestry( adjacency )
  
  -sum( getScoreBounds( lll, ancestors, !ancestors )["upper",] )
}

#' Toggle search children
#' 
#' Powers an edge-toggling search over adjacency matrices
#' @param adjacency
#' @return list of adjacency matrices resulting from toggling an edge
#' @export
getChildrenInToggleSearch = function ( adjacency ){
  n = length( adjacency )
  
  lapply( 1:n, function( x ) {
    xor( adjacency, logicVector( n, x ) )
  } )
}