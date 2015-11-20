# search-with-splits.R
# Searches for a best network considering splits.

# This version looks for the best actor ancestry individually per reporter and then
# composes.
searchWithSplits = function( laps, searchFunction = searchByAStar){
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

# This version is a variation on A* search
searchByAStar = function( laps, reporterIndex ){
  
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
    index = which.max(boundary.scores)
    score = boundary.scores[ index ]
    state = boundary[[ index ]]
    boundary = boundary[ -index ]
    boundary.scores = boundary.scores[ -index ]
    
    if ( state$depth == maxDepth ){
      ancestry = deriveAncestry(
        state$adjacency,
        matrix( FALSE, nrow = n, ncol = n + 1 ) )$ancestry
      break
    }
    
    # We can reuse index as our insertion point now
    index = length( boundary ) + 1
    
    # Figure out our next edge toggle
    toggle.point = which( cbind( diag( n ) == 0, TRUE ) )[ n^2 - state$depth ]
    
    # Compute "edge present"
    new.adjacency = state$adjacency
    new.adjacency[ toggle.point ] = TRUE
    new.state = list( depth = state$depth + 1, adjacency = new.adjacency )
    # "push" new state in.
    boundary[ index ] = list( new.state )
    boundary.scores[ index ] = getHeuristicScore(
      laps,
      reporterIndex,
      state$depth + 1,
      new.adjacency )
    
    index = index + 1
    # Compute "edge absent"
    new.adjacency[ toggle.point ] = FALSE
    new.state = list( depth = state$depth + 1, adjacency = new.adjacency )
    # "push" new state in.
    boundary[ index ] = list( new.state )
    boundary.scores[ index ] = getHeuristicScore(
      laps,
      reporterIndex,
      state$depth + 1,
      new.adjacency )
  }
  
  return ( list( score = score, ancestry = ancestry ) )
}

# The heuristic + score for A*
getHeuristicScore = function ( laps, reporterIndex, depth, adjacency ){
  
  n = nrow( adjacency )
  
  # Mark certain/free edges in adjacency matrix.
  uncertain = cbind( diag( n ) == 0, TRUE )
  
  if( depth > 0 ){
    indeces = which( uncertain )
    uncertain[ indeces[ length( indeces ) - ( ( depth - 1 ):0 ) ] ] = FALSE
  }
  
  #print( uncertain )
  
  # Derive ancestry
  ancestral = deriveAncestry( adjacency, uncertain )
  uncertain = ancestral$uncertain
  ancestry = ancestral$ancestry
  
  #print( uncertain )
  #print( ancestry )
  
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
  
  return ( score )
}

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