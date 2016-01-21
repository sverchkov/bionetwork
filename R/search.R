# Score search

#' Search for a best network (without actor duplication)
#' 
#' This version is a variation on A* search
#' 
#' @param lll the LocalLogLikelihoods object
#' @param reporters the reporters for which we're building the network (all reporters in lll by default)
#' @param searchSelection function determining the search order, by default using lastMax,
#' which corresponds to DFS A*. An interesting alternative is choosing a random max,
#' choosing the first max will yield the abysmally slow BFS A*.
#' @param debug.level granularity of debug statements. 0 is none, larger numbers mean more statements.
#' @return a list containing the ancestry structure specification (two adjacency matrices: one for present edges, and one for uncertain edges) and the reporter-wise score vector.
#' @export
branchBoundAStarOnTree = function( lll, reporters = getReporters( lll ), searchSelection = whichLastMax, debug.level = 0 ){
  
  n = howManyActors( lll ) # Number of actors
  nR = length( reporters ) # Number of reporters
  
  the.dim.names = list( getActors( lll ), c( getActors( lll ), reporters ) )
  
  returned.score = c( -Inf, -Inf )
  adjacency = matrix( FALSE, nrow = n, ncol = n + nR, dimnames = the.dim.names )
  uncertain = ( 1 != diag( nrow = n, ncol = n + nR ) )
  dimnames( uncertain ) = the.dim.names
  best.lowerbound = -Inf
  
  # Initial search insertion
  boundary = list( list(
    adjacency = adjacency,
    uncertain = uncertain
  ) ) # Start with a single dumb-state
  scores = matrix( c( Inf, -Inf ), nrow = 2, ncol = 1, dimnames = list( c("upper","lower") ) )
  
  repeat {
    
    if( debug.level > 5 ){
      print( "Boundary size:" )
      print( length( boundary ) )
      print( "Score matrix size:" )
      print( dim( scores ) )
    }
    
    # Find and "pop" best scoring point
    index = searchSelection( scores["upper",] )
    
    if( debug.level > 5 ) print( index )
    
    score = scores[, index ]
    state = boundary[[ index ]]
    boundary = boundary[ -index ]
    scores = scores[, -index ]
    
    uncertain = state$uncertain
    adjacency = state$adjacency

    # Stop condition: When upperbound = lowerbound we can't make any useful decisions anymore
    if ( score["lower"] >= score["upper"] ){
      
      # Build uncertainty-based ancestry matrices
      possible.ancestors = adjacencyToAncestry( new.adjacency | new.uncertain )
      possible.nonancestors = !adjacencyToAncestry( new.adjacency & !new.uncertain )
      
      # Get that score vector for return
      returned.score = getScoreBounds( lll, possible.ancestors, possible.nonancestors )
      
      break
    }
    
    # Ok, we're going forward, so let's update our best lowerbound
    best.lowerbound = max( best.lowerbound, score["lower"] )
    
    for ( new.state in getNextSearchStates( adjacency, uncertain, last.edge = FALSE ) ){

      if( debug.level > 4 ) print( "Making branch..." )
      
      # Build uncertainty-based ancestry matrices
      possible.ancestors = adjacencyToAncestry( new.state$adjacency | new.state$uncertain )
      possible.nonancestors = !adjacencyToAncestry( new.state$adjacency & !new.state$uncertain )
      
      # Score
      new.score = rowSums( getScoreBounds( lll, possible.ancestors, possible.nonancestors ) )

      if( debug.level > 3 ) print( new.score )
      
      # Bound check pruning
      if( new.score["upper"] >= best.lowerbound ) {
        # "push" new state in.
        boundary[[ length(boundary) + 1 ]] = new.state
        scores = cbind( scores, new.score )
        if( debug.level > 4 ) print( "Added branch." )
      }
    }
  }
  
  # For debugging
  if ( debug.level > 2 ) print( returned.score )
  
  return ( list( score = returned.score, adjacency = adjacency, uncertain = uncertain ) )
}

#' Search for a best network (without actor duplication)
#' 
#' This version is a variation on A* search
#' 
#' @param lll the LocalLogLikelihoods object
#' @param reporters the reporters for which we're building the network (all reporters in lll by default)
#' @param searchSelection function determining the search order, by default using lastMax,
#' which corresponds to DFS A*. An interesting alternative is choosing a random max,
#' choosing the first max will yield the abysmally slow BFS A*.
#' @param debug.level granularity of debug statements. 0 is none, larger numbers mean more statements.
#' @return a list containing the ancestry structure specification (two adjacency matrices: one for present edges, and one for uncertain edges) and the reporter-wise score vector.
#' @export
branchBoundAStarOnGraph = function( lll, reporters = getReporters( lll ), searchSelection = whichLastMax, debug.level = 0 ){
  
  n = howManyActors( lll ) # Number of actors
  nR = length( reporters ) # Number of reporters
  
  the.dim.names = list( getActors( lll ), c( getActors( lll ), reporters ) )
  
  returned.score = c( -Inf, -Inf )
  adjacency = matrix( FALSE, nrow = n, ncol = n + nR, dimnames = the.dim.names )
  uncertain = ( 1 != diag( nrow = n, ncol = n + nR ) )
  dimnames( uncertain ) = the.dim.names
  best.lowerbound = -Inf
  
  # Initial search insertion
  boundary = list( list(
    adjacency = adjacency,
    uncertain = uncertain
  ) ) # Start with a single dumb-state
  scores = matrix( c( Inf, -Inf ), nrow = 2, ncol = 1, dimnames = list( c("upper","lower") ) )
  
  repeat {
    
    if( debug.level > 5 ){
      print( "Boundary size:" )
      print( length( boundary ) )
      print( "Score matrix size:" )
      print( dim( scores ) )
    }
    
    # Find and "pop" best scoring point
    index = searchSelection( scores["upper",] )
    
    if( debug.level > 5 ) print( index )
    
    score = scores[, index ]
    state = boundary[[ index ]]
    boundary = boundary[ -index ]
    scores = scores[, -index ]
    
    uncertain = state$uncertain
    adjacency = state$adjacency
    
    # Stop condition: When upperbound = lowerbound we can't make any useful decisions anymore
    if ( score["lower"] >= score["upper"] ){
      
      # Build uncertainty-based ancestry matrices
      possible.ancestors = adjacencyToAncestry( new.adjacency | new.uncertain )
      possible.nonancestors = !adjacencyToAncestry( new.adjacency & !new.uncertain )
      
      # Get that score vector for return
      returned.score = getScoreBounds( lll, possible.ancestors, possible.nonancestors )
      
      break
    }
    
    # Ok, we're going forward, so let's update our best lowerbound
    best.lowerbound = max( best.lowerbound, score["lower"] )
    
    for ( new.state in getNextGraphSearchStates( adjacency, uncertain ) ){
      
      if( debug.level > 4 ) print( "Making branch..." )
      
      # Build uncertainty-based ancestry matrices
      possible.ancestors = adjacencyToAncestry( new.state$adjacency | new.state$uncertain )
      possible.nonancestors = !adjacencyToAncestry( new.state$adjacency & !new.state$uncertain )
      
      # Score
      new.score = rowSums( getScoreBounds( lll, possible.ancestors, possible.nonancestors ) )
      
      if( debug.level > 3 ) print( new.score )
      
      # Bound check pruning
      if( new.score["upper"] >= best.lowerbound ) {
        # "push" new state in.
        boundary[[ length(boundary) + 1 ]] = new.state
        scores = cbind( scores, new.score )
        if( debug.level > 4 ) print( "Added branch." )
      }
    }
  }
  
  # For debugging
  if ( debug.level > 2 ) print( returned.score )
  
  return ( list( score = returned.score, adjacency = adjacency, uncertain = uncertain ) )
}

#' Make uncertainty matrix from depth + dimension
#' 
#' @param depth - search depth
#' @param n - number of rows
#' @param m - number of columns
getUncertaintyMatrix = function ( depth = 0, n, m = n ){
  
  # Diagonals are never uncertain.
  uncertain = ( diag( nrow = n, ncol = m ) == 0 )
  
  if( depth > 0 ){
    indeces = which( uncertain )
    uncertain[ indeces[ length( indeces ) - ( ( depth - 1 ):0 ) ] ] = FALSE
  }
  
  return ( uncertain )
}

#' Get the next state branches in the search
#' 
#' If the possible ancestry/nonancestry matrices are provided, the process is
#' sped up in some cases.
#' 
#' @param adjacency the present adjacency matrix
#' @param uncertain the present uncertain matrix
#' @param last.edge whether to pick the last uncertain edge to make certain (default is
#' true), if false picks the first uncertain edge instead.
#' @return a list of two elements (corresponding to making the last (or first)
#' ucertain edge certain) where each is the state description of the next branch,
#' consisting of an adjacency and an uncertain edge matrix
getNextSearchStates = function( adjacency
                              , uncertain
                              , last.edge = TRUE ){
  if ( !any( uncertain ) )
    stop("Something is wrong: trying to search, but no edges are uncertain!")
  
  n = nrow( adjacency )
  
  ind = which( uncertain, arr.ind = TRUE )
  i = ind[ nrow( ind ), 1 ]
  j = ind[ nrow( ind ), 2 ]
  if ( !last.edge ){
    i = ind[ 1, 1 ]
    j = ind[ 1, 2 ]
  }
  
  new.uncertain = uncertain
  new.uncertain[ i, j ] = FALSE
  
  result = list()
  
  for ( new.edge in c( FALSE, TRUE ) ){
    new.adjacency = adjacency
    new.adjacency[ i, j ] = new.edge
    
    result = c( result, list( list( adjacency = new.adjacency, uncertain = new.uncertain ) ) )
  }
  
  result
}

#' Get the next state branches a graph search
#' 
#' @param adjacency the present adjacency matrix
#' @param uncertain the present uncertain matrix
#' @param fragile whether to fail when there are no uncertain edges, (returns an empty
#' list if not fragile)
#' @return a list of two elements (corresponding to making the last (or first)
#' ucertain edge certain) where each is the state description of the next branch,
#' consisting of an adjacency and an uncertain edge matrix
getNextGraphSearchStates = function( adjacency
                                   , uncertain
                                   , fragile = TRUE ){
  if ( fragile && !any( uncertain ) )
    stop("Something is wrong: trying to search, but no edges are uncertain!")
  
  result = list()
  
  for ( i in which( uncertain ) ) {

    new.uncertain = uncertain
    new.uncertain[ i ] = FALSE
  
    for ( new.edge in c( FALSE, TRUE ) ){
      new.adjacency = adjacency
      new.adjacency[ i ] = new.edge
      
      result = c( result, list( list( adjacency = new.adjacency, uncertain = new.uncertain ) ) )
    }
  }
  
  result
}