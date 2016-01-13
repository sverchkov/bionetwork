# Score functions

#' Likelihood scores broken down by reporter
#' 
#' Likelihood calculation for a full network, with the score contribution of each
#' reporter returned in a vector.
#' @param ancestry - the ancestry matrix
#' @param lll - the LocalLogLikelihoods object
#' @return a vector with a log-likelihood for each reporter
scoreLikelihoodsPerReporter = function ( ancestry, lll ) {
  
  actor.list = rownames( ancestry )
  nA = length( actor.list )
  
  nR = howManyReporters( lll )
  unique.actors = getActors( lll )
  
  # For a reporter r
  # Ancestors of r :  ancestry[, r ]
  # Named :           actor.list[ ancestry[, r ] ]
  # Named as in lll : stripDots( actor.list[ ancestry[, r ] ] )
  # Selected in lll : unique.actors %in% stripDots( actor.list ... )
  non.ancestors = sapply( as.list( getReporters( lll ) ), function( r ) {
    unique.actors %in% stripDots( actor.list[ !ancestry[, r ] ] )
  } )
  
  # Simple ancestry component:
  simple.ancestry.scores = ancestryScoreMatrix( lll )
  simple.ancestry.scores[ non.ancestors ] = nonAncestryScoreMatrix( lll )[ non.ancestors ]
  
  # Clean NAs
  simple.ancestry.scores[ is.na( simple.ancestry.scores ) ] = 0
  
  # Get scores
  scores = colSums( simple.ancestry.scores )
  
  # Scoring 2le KOs requires some special care. For "neither ancestor" we need to make
  # sure we aren't scoring a reporter multiple times, so we need to iterate over the
  # unique actors. For the other states we use the pooled versions (and as long as our
  # starting ancestry is correct there shouldn't be any double-counting)
  
  # 2le KO "neither ancestor" component:
  for ( a in 2:length( unique.actors ) ){
    for ( b in 1:a ){
      selected.reporters = !apply(
        ancestry[ stripDots( actor.list ) %in% unique.actors[ c(a,b) ],
                  nA + ( 1:nR ) ],
        2, any )
      scores[ selected.reporters ] = scores[ selected.reporters ] + replaceNAs( scoreNeitherAncestor( lll, unique.actors[ a ], unique.actors[ b ] )[ selected.reporters ] )
    }
  }
  
  # 2le KO components for remainder:
  for ( a in 2:nA ){
    actor.a.pooled = actor.list[a]
    actor.a.unique = stripDots( actor.a.pooled )
    for ( b in 1:a ){
      actor.b.pooled = actor.list[b]
      actor.b.unique = stripDots( actor.b.pooled )

      # Ancestor masks:
      # ancestry[ a, nA + ( 1:nR ) ] <- ancestors of a
      only.a = ancestry[ a, nA + (1:nR) ] & !ancestry[ b, nA + (1:nR) ]
      only.b = ancestry[ b, nA + (1:nR) ] & !ancestry[ a, nA + (1:nR) ]
      both = ancestry[ a, nA + (1:nR) ] & ancestry[ b, nA + (1:nR) ]
      
      # Only A ancestor score
      scores[ only.a ] = scores[ only.a ] + replaceNAs( scoreSingleAncestor( lll, actor.a.unique, actor.b.unique )[ only.a ] )
      # Only B ancestor score
      scores[ only.b ] = scores[ only.b ] + replaceNAs( scoreSingleAncestor( lll, actor.b.unique, actor.a.unique )[ only.b ] )
      # Both
      scores[ both ] = scores[ both ] + replaceNAs(
        if ( ancestry[ a, b ] || ancestry[ b, a ] ){
          # Shared pathway score
          scoreSharedPathways( lll, actor.a.unique, actor.b.unique )[ both ]
        } else {
          # Independent pathway score
          scoreIndependentPathways( lll, actor.a.unique, actor.b.unique )[ both ]
        } )
    }
  }
      
  return ( scores )
}

#' The heuristic + score for A*
#' This one uses the log-likelihood interface
getHeuristicScore = function ( lll, reporterIndex, depth, adjacency ){
  
  n = nrow( adjacency )
  actors = getActors( lll )
  
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
    sum( pmax(
      ancestryScoreMatrix( lll )[ uncertain[, n + 1 ], reporterIndex ],
      nonAncestryScoreMatrix( lll )[ uncertain[, n + 1 ], reporterIndex ] ), na.rm = TRUE )
  
  # Sorta hacky
  if ( is.na( score ) ) score = 0
  
  # 2le KO component:
  for( a in 2:n ){
    for ( b in 1:a ){
      rel.score = c(
        # Neither ancestor score
        neither = scoreNeitherAncestor( lll, actors[a], actors[b] )[ reporterIndex ],
        # Only A ancestor score
        only.a = scoreSingleAncestor( lll, actors[a], actors[b] )[ reporterIndex ],
        # Only B ancestor score
        only.b = scoreSingleAncestor( lll, actors[b], actors[a] )[ reporterIndex ],
        # Both, independent pathway score
        independent = scoreIndependentPathways( lll, actors[a], actors[b] )[ reporterIndex ],
        # Both, shared pathway score
        shared = scoreSharedPathways( lll, actors[a], actors[b] )[ reporterIndex ] )
      
      # Sorta hacky
      rel.score[ is.na( rel.score ) ] = 0
      
      # Now go through the cases
      if ( !uncertain[ a, n + 1 ] ){ # Only certain relations eliminate cases
        if ( ancestry[ a, n + 1 ] ) { # a is ancestor
          rel.score[ c( "only.b", "neither" ) ] = -Inf
        } else {
          rel.score[ c( "only.a", "independent", "shared" ) ] = -Inf
        }
      }
      if ( !uncertain[ b, n + 1 ] ){
        if ( ancestry[ b, n + 1 ] ){ # b is ancestor
          rel.score[ c( "only.a", "neither" ) ] = -Inf
        } else { # b is not ancestor
          rel.score[ c( "only.b", "independent", "shared" ) ] = -Inf
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
