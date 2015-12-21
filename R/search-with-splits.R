# search-with-splits.R
# Searches for a best network considering splits.

#' Search with splits
#' 
#' This version looks for the best actor ancestry individually per reporter and then
#' composes.
searchWithSplits = function( laps, searchFunction = searchByAStar ){
  ancestry.template = matrix(
    FALSE,
    nrow = howManyActors( laps ),
    ncol = howManyActors( laps ) + 1,
    dimnames = list( getActors( laps ), c( getActors( laps ), "r" ) ) )
  
  results.list = lapply(
    as.list( 1:howManyReporters( laps ) ),
    function( x ) searchFunction( laps, x ) )
  
  scores = vapply( results.list, function ( x ) x$score, -Inf )
  ancestries = vapply( results.list, function ( x ) x$ancestry, ancestry.template )
    
  return ( list( scores = scores, ancestries = ancestries ) )
}

#' Search for a reporter's best network
#' 
#' This version is a variation on A* search
#' 
#' @param laps the Local Ancestry Probablity Score object
#' @param reporterIndex the index of the reporter for which we're building the network
#' @param searchSelection function determining the search order, by default using lastMax,
#' which corresponds to DFS A*. An interesting alternative is choosing a random max,
#' choosing the first max will yield the abysmally slow BFS A*.
searchByAStar = function( laps, reporterIndex, searchSelection = whichLastMax ){
  
  # For debugging
  print( reporterIndex )
  
  n = howManyActors( laps )
  
  maxDepth = n^2
  
  score = -Inf
  adjacency = matrix( FALSE, nrow = n, ncol = n + 1 )
  
  # Initial search insertion
  boundary = list( list(
    depth = 0,
    adjacency = matrix( FALSE, nrow = n, ncol = n + 1 )
  ) ) # Start with a single dumb-state
  boundary.scores = 0
  
  repeat {
    # Find and "pop" best scoring point
    index = searchSelection( boundary.scores )
    score = boundary.scores[ index ]
    state = boundary[[ index ]]
    boundary = boundary[ -index ]
    boundary.scores = boundary.scores[ -index ]
    
    # For debugging
    #print('Boundary:')
    #print(boundary.scores)
    #print('Depth:')
    #print(state$depth)
    
    uncertainAncestry = deriveAncestry(
      state$adjacency,
      getUncertaintyMatrix( state$depth, n, n+1 ) )
    
    if ( state$depth == maxDepth || !any( uncertainAncestry$uncertain[,n+1] ) ){
      ancestry = uncertainAncestry$ancestry
      break
    }
    
    # For debugging
    #print('Matrix:')
    #print( showUncertainAncestry( uncertainAncestry$ancestry, uncertainAncestry$uncertain ) )
    
    # We can reuse index as our insertion point now
    index = length( boundary ) + 1
    
    # Figure out our next edge toggle
    toggle.point = which( cbind( diag( n ) == 0, TRUE ) )[ n^2 - state$depth ]
    
    # Compute "edge present"
    new.adjacency = state$adjacency
    new.adjacency[ toggle.point ] = TRUE
    new.state = list( depth = state$depth + 1, adjacency = new.adjacency )
    # pruning and scoring
    if( -Inf < ( new.score = getHeuristicScore( laps,
                                                reporterIndex,
                                                state$depth + 1,
                                                new.adjacency ) ) ) {
      # "push" new state in.
      boundary[ index ] = list( new.state )
      boundary.scores[ index ] = new.score
      index = index + 1
    }
    
    # Compute "edge absent"
    new.adjacency[ toggle.point ] = FALSE
    new.state = list( depth = state$depth + 1, adjacency = new.adjacency )
    # "push" new state in.
    if( -Inf < ( new.score = getHeuristicScore( laps,
                                                reporterIndex,
                                                state$depth + 1,
                                                new.adjacency ) ) ) {
      # "push" new state in.
      boundary[ index ] = list( new.state )
      boundary.scores[ index ] = new.score
    }
  }
  
  # For debugging
  print( score )
  
  return ( list( score = score, ancestry = ancestry ) )
}

#' Given a vector, get the index of the last maximal element
#' 
#' @param a vector of numbers
#' @return the index of the last maximal element
whichLastMax = function ( x ){
  z = which( x == max( x ) )
  z[length(z)]
}

#' Given a vector, get the index of a random maximal element
#' 
#' Given a vector, looks at its maximal elements, walking from last to first,
#' rejecting each with probability (1-p). (So has probability p to pick the last,
#' p*(1-p)th the 2nd last, etc. with the remainder of the probabilty mass placed
#' in the first).
#' 
#' @param a vector of numbers
#' @return the index of a (non-uniformly) random maximal element
whichRandomLateMax = function ( x, p = 0.9 ){
  z = which( x == max( x ) )
  i = length(z)
  while( sample( x = c( TRUE, FALSE ), size = 1, prob = c( 1-p, p ) ) && i > 1 ){
    i = i-1
  }
  z[i]
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

#' The heuristic + score for A*
#' This one uses the log-likelihood interface
getHeuristicScore = function ( lll, reporterIndex, depth, adjacency ){
  
  n = nrow( adjacency )
  acotrs = getActors( lll )
  
  # Mark certain/free edges in adjacency matrix.
  uncertain = getUncertaintyMatrix( depth, n, n+1 )
  
  #print( uncertain )
  
  # Derive ancestry
  ancestral = deriveAncestry( adjacency, uncertain )
  uncertain = ancestral$uncertain
  ancestry = ancestral$ancestry
  
  # For debugging  
  #print( showUncertainAncestry( ancestry, uncertain ) )
  
  # Score
  
  # Simple ancestry component.
  score =
    # Certain ancestors
    sum( ancestryScoreMatrix( lll )[ ancestry[, n + 1 ] & !uncertain[, n + 1 ], reporterIndex], na.rm = TRUE ) +
    # Certain non-ancestors
    sum( nonAncestryScoreMatrix( lll )[ !ancestry[, n + 1 ] & !uncertain[, n + 1 ], reporterIndex ], na.rm = TRUE ) +
    # Uncertain
    sum( mapply(
      max,
      ancestryScoreMatrix( lll )[ uncertain[, n + 1 ], reporterIndex ],
      nonAncestryScoreMatrix( lll )[ uncertain[, n + 1 ], reporterIndex ] ), na.rm = TRUE )
  
  # 2le KO component:
  for( a in 2:n ){
    for ( b in 1:a ){
      rel.score = c(
        # Neither ancestor score
        neither = scoreNeitherAncestor( lll, actors(a), actors(b) )[ reporterIndex ],
        # Only A ancestor score
        only.a = scoreSingleAncestor( lll, actors(a), actors(b) )[ reporterIndex ],
        # Only B ancestor score
        only.b = scoreSingleAncestor( lll, actors(b), actors(a) )[ reporterIndex ],
        # Both, independent pathway score
        independent = scoreIndependentPathways( lll, actors(a), actors(b) )[ reporterIndex ],
        # Both, shared pathway score
        shared = scoreSharedPathways( lll, actors(a), actors(b) )[ reporterIndex ] )
      
      # Now go through the cases
      if ( !uncertain[ a, n + 1 ] ){ # Only certain relations eliminate cases
        if ( ancestry[ a, n + 1 ] ) { # a is ancestor
          rel.score[ "only.b", "neither" ] = -Inf
        } else {
          rel.score[ "only.a", "independent", "shared" ] = -Inf
        }
      }
      if ( !uncertain[ b, n + 1 ] ){
        if ( ancestry[ b, n + 1 ] ){ # b is ancestor
          rel.score[ "only.a", "neither" ] = -Inf
        } else { # b is not ancestor
          rel.score[ "only.b", "independent", "shared" ] = -Inf
        }
      }
      if ( ( !uncertain[ a, b ] && ancestry[ a, b ] ) || ( !uncertain[ b, a ] && ancestry[ b, a ] ) ) {
        rel.score[ "independent" ] = -Inf
      } else if ( !uncertain[ a, b ] && !uncertain[ b, a ] && !ancestry[ a, b ] && !ancestry[ b, a ] ) {
        rel.score[ "shared" ] = -Inf
      }
      
      score = score + max( rel.score )
    }
  }
  
  # For debugging
  #print( score )
  
  return ( score )
}

#' The heuristic + score for A*
getHeuristicScoreOld = function ( laps, reporterIndex, depth, adjacency ){
  
  n = nrow( adjacency )
  
  # Mark certain/free edges in adjacency matrix.
  uncertain = getUncertaintyMatrix( depth, n, n+1 )
  
  #print( uncertain )
  
  # Derive ancestry
  ancestral = deriveAncestry( adjacency, uncertain )
  uncertain = ancestral$uncertain
  ancestry = ancestral$ancestry

  # For debugging  
  #print( showUncertainAncestry( ancestry, uncertain ) )
  
  # Score
  
  # Simple ancestry component:
  local.scores = ancestryScoreMatrix( laps )[,reporterIndex]
  score =
    # Certain ancestors
    sum( local.scores[ ancestry[, n + 1 ] & !uncertain[, n + 1 ] ], na.rm=TRUE ) +
    # Certain non-ancestors
    sum( log1mexp(
      local.scores[ !ancestry[, n + 1 ] & !uncertain[, n + 1 ] ] ), na.rm=TRUE ) +
    # Uncertain
    if ( any ( uncertain[, n + 1 ] ) ){
      sum( mapply(
        max,
        local.scores[ uncertain[, n + 1 ] ],
        log1mexp( local.scores[ uncertain[, n + 1 ] ] ) ), na.rm = TRUE )
    }else 0
  
  # 2le KO component:
  actors = which( ancestry[, n + 1 ] | uncertain[, n + 1 ] )
  if ( length( actors ) > 1 ){
    for ( i in 2:length( actors ) ){
      for ( b in actors[ 1:( i - 1 ) ] ){
        a = actors[ i ]
        sp.score = scoreSharedPathways( laps, a, b )[reporterIndex]
        ip.score = scoreIndependentPathways( laps, a, b )[reporterIndex]
        #print( score )
        #print( sp.score )
        #print( ip.score )
        score = score +
          if ( any( uncertain[ c( a, b ), c( a, b, n + 1 ) ] ) ){
            max( sp.score, ip.score )
          } else {
            if ( ancestry[a,b] || ancestry[b,a] )
              sp.score
            else
              ip.score
          }
      }
    }
  }
  
  # For debugging
  #print( score )
  
  return ( score )
}

#' Derive ancestry from uncertain adjacency
deriveAncestry = function ( adjacency, uncertain ){
  # Derives an ancestry from a partially uncertain adjacency matrix
  
  # Get our size (assume adjacency dims = uncertain dims )
  n = nrow( adjacency )
  m = ncol( adjacency )
  
  # First, ensure that only certain edges are in the adj matrix
  adjacency[ uncertain ] = FALSE
  
  # Next, derive certain ancestry
  ancestry = ( diag( m ) == TRUE )
  updated.ancestry = matrix( FALSE, nrow = m, ncol = m )
  updated.ancestry[1:n,1:m] = ancestry[1:n,1:m] | adjacency
  while ( any( ancestry != updated.ancestry ) ){
    ancestry = updated.ancestry
    for ( node in 1:m ){
      updated.ancestry[, node ] =
        apply( as.matrix( ancestry[, ancestry[, node ] ] ), 1, any )
    }
  }
  
  updated.uncertain = matrix( FALSE, nrow = m, ncol = m )
  updated.uncertain[1:n,1:m] = uncertain
  uncertain = matrix( FALSE, nrow = m, ncol = m )
  while ( any( uncertain != updated.uncertain ) ){
    uncertain = updated.uncertain
    propagator = updated.uncertain | updated.ancestry
    for ( node in 1:m ){
      updated.uncertain[, node ] =
        apply( as.matrix( updated.uncertain[, propagator[, node ] ] ), 1, any )
    }
  }
  updated.uncertain[updated.ancestry] = FALSE
  
  return ( list( ancestry = updated.ancestry[1:n,1:m], uncertain = updated.uncertain[1:n,1:m] ) )
}

#' A human-friendly representation of ancestry+uncertainty matrix
#' 
#' @param ancestry the ancesty matrix
#' @param uncertainty the uncertainty matrix
#' @return a char matrix with T/F in certain cells and ? in uncertain cells.
showUncertainAncestry = function ( ancestry, uncertainty ){
  dims = dim( ancestry )
  result = matrix( '?', nrow = dims[1], ncol = dims[2] )
  result[ ancestry & !uncertainty ] = 'T'
  result[ !ancestry & !uncertainty ] = 'F'
  return ( result )
}