# A method for inferring the maximum-likelihood network based on log-odds from
# contrasts.
#
# Initial creation: August 10, 2015
# Yuriy Sverchkov

# Interface used by the network learner.
# I'm guessing in R the interface is just generics.
setGeneric(name = "scoreIndependentPathways",
           def = function(theObject, actor1, actor2)
             { standardGeneric("scoreIndependentPathways") } )
setGeneric(name = "scoreSharedPathways",
           def = function(theObject, actor1, actor2)
             { standardGeneric("scoreSharedPathways") } )
setGeneric(name = "ancestryScoreMatrix",
           def = function(theObject)
             { standardGeneric("ancestryScoreMatrix") } )
setGeneric(name = "getActors",
           def = function(theObject)
           { standardGeneric("getActors") } )
setGeneric(name = "getReporters",
           def = function(theObject)
           { standardGeneric("getReporters") } )
setGeneric(name = "howManyActors",
           def = function(theObject)
           { standardGeneric("howManyActors") } )
setGeneric(name = "howManyReporters",
           def = function(theObject)
             { standardGeneric("howManyReporters") } )

# The class LimmaLogProbs is defined in the "using-limma.R" source file.
# Here we associate the generics with it.
setMethod(f = "scoreIndependentPathways",
          signature(theObject = "LimmaLogProbs", actor1 = "character", actor2 = "character"),
          definition = function(theObject, actor1, actor2){
            (
            if( actor1 %in% names(theObject@doubleVsingle) &&
                actor2 %in% names( theObject@doubleVsingle[[actor1]] ) )
                theObject@doubleVsingle[[actor1]][[actor2]][["gt"]]
            else -log(2)
            ) + (
            if( actor2 %in% names(theObject@doubleVsingle) &&
                actor1 %in% names( theObject@doubleVsingle[[actor2]] ) )
                theObject@doubleVsingle[[actor2]][[actor1]][["gt"]]
            else -log(2)
            )
          })

setMethod(f = "scoreIndependentPathways",
          signature(theObject = "LimmaLogProbs", actor1 = "numeric", actor2 = "numeric"),
          definition = function(theObject,actor1,actor2){
            scoreIndependentPathways(theObject, theObject@actors[actor1],theObject@actors[actor2])
          })

setMethod(f = "scoreSharedPathways",
          signature(theObject = "LimmaLogProbs", actor1 = "character", actor2 = "character"),
          definition = function(theObject, actor1, actor2){
            (
            if( actor1 %in% names(theObject@doubleVsingle) &&
                actor2 %in% names( theObject@doubleVsingle[[actor1]] ) )
                theObject@doubleVsingle[[actor1]][[actor2]][["eq"]]
            else -log(2)
            ) + (
            if( actor2 %in% names(theObject@doubleVsingle) &&
                actor1 %in% names( theObject@doubleVsingle[[actor2]] ) )
                theObject@doubleVsingle[[actor2]][[actor1]][["eq"]]
            else -log(2)
            )
          })

setMethod(f = "scoreSharedPathways",
          signature(theObject = "LimmaLogProbs", actor1 = "numeric", actor2 = "numeric"),
          definition = function(theObject,actor1,actor2){
            scoreSharedPathways(theObject, theObject@actors[actor1],theObject@actors[actor2])
          })

setMethod(f = "ancestryScoreMatrix",
          signature = "LimmaLogProbs",
          definition = function(theObject) t(theObject@singleGtWT) )

setMethod(f = "getActors",
          signature = "LimmaLogProbs",
          definition = function(theObject) theObject@actors )

setMethod(f = "getReporters",
          signature = "LimmaLogProbs",
          definition = function(theObject) theObject@reporters )

setMethod(f= "howManyActors",
          signature = "LimmaLogProbs",
          definition = function(theObject) theObject@nActors )

setMethod(f= "howManyReporters",
          signature = "LimmaLogProbs",
          definition = function(theObject) theObject@nReporters )

# Multiple start (by default greedy) search,
# only over ancesty-like networks:
multiStartANetworkSearch <- function(
  lp,
  searchFunction=getMLNetwork)
  # See getMLNetwork for parameter definitions
{
  # Let n = number of actors
  actors = getActors(lp)
  n = howManyActors(lp)
  nReporters = howManyReporters(lp)
  # We want to do a search starting with each combination of the possible
  # n^2 - n edges.
  # Start by getting the list of edges
  edges = which( matrix( TRUE, nrow=n, ncol=n ) & ! diag(n) )
  
  # Keep track of max score
  max.score = -Inf
  max.network = NULL
  
  # Next iterate over 0-max num of edges
  for( n.edges in 0:length(edges) )
    # Iterate over combinations
    for( combination in combn(edges, n.edges, simplify = FALSE) ){
      # Make ancestry graph
      ancestry = diag(n)
      for( i in combination ) ancestry[i] = TRUE
      
      # Check ancestry (skip if not holding)
      if ( all( ancestry == inclusive.ancestry( ancestry ) ) ){
  
        # Get a network from the ancestry graph      
        network = cbind( ancestry & !diag(n), matrix( FALSE, nrow=n, ncol=nReporters ) )
        
        ml.result =
          searchFunction( lp, n, nReporters, network )
        
        candidate.network = ml.result[[1]]
        candidate.score = ml.result[[2]]
        if( candidate.score > max.score ){
          max.score = candidate.score
          max.network = candidate.network
        }
      }
    }
  
  return(list(max.score,max.network))
}

# Multiple start (by default greedy) search:
multiStartNetworkSearch <- function(
  lp,
  searchFunction=getMLNetwork)
  # See getMLNetwork for parameter definitions
{
  # Let n = number of actors
  actors = getActors(lp)
  n = howManyActors(lp)
  nReporters = howManyReporters(lp)
  # We want to do a search starting with each combination of the possible
  # n^2 - n edges.
  # Start by getting the list of edges
  edges = which(
    matrix( TRUE, nrow=n, ncol=n ) &
    !diag( n ) )
  
  # Keep track of max score
  max.score = -Inf
  max.network = NULL
  
  # Next iterate over 0-max num of edges
  for( n.edges in 0:length(edges) )
    # Iterate over combinations
    for( combination in combn(edges, n.edges, simplify = FALSE) ){
      # Make a starting network
      network = matrix( FALSE, nrow=n, ncol=(n+nReporters) )
      for( i in combination ) network[i] = TRUE
      
      # Score
      #candidate.score = scoreNetwork( network, lp )
      
      ml.result =
        searchFunction( lp, n, nReporters, network )
    
      candidate.network = ml.result[[1]]
      candidate.score = ml.result[[2]]
      if( candidate.score > max.score ){
        max.score = candidate.score
        max.network = candidate.network
      }
    }
  
  return(list(max.score,max.network))
}

# Greedy edge toggle search:
getMLNetwork <- function(
  lp,
  nActors = length(getActors(lp)),
  nReporters = howManyReporters(lp),
  network=matrix( FALSE, nrow=nActors,  ncol=(nActors+nReporters)),
  score=scoreNetwork( network, lp ) )
  # lods: data structure describing the log-odds comparing various KOs
  # n.actors: the number of things upstream in the network
  # n.reporters: the number of genes affected by the actors
  # network: initial network from which to begin the search (default is empty)
  # score: initial score of initial network
  # Graph is represented by the directed adjacency matrix where network[a,b] == TRUE
  # indicates the edge a->b
{
  print( paste0( nActors, " actors and ", nReporters, " reporters." ))
  repeat{
    print(paste0("Score: ",score))
    score.change.matrix <- scoreEdgeToggles( network, lp )
    max.change <- max(score.change.matrix, na.rm = TRUE)
    if ( max.change <= 0 || max.change == Inf ) {
      return(list(network,score,score.change.matrix))
    } else {
      while ( max.change > 0 ) {
        max.index <- which(score.change.matrix == max.change, arr.ind = TRUE)
        row = sample( 1:(dim(max.index)[1]), 1 )
        a = max.index[row,1]
        b = max.index[row,2]
        print( paste0( "Toggled ", a, ",", b, " yielding change: ", max.change ))
        network[a,b] = !network[a,b]
        score = score + max.change
        if ( b <= nActors ) break
        # This is a short-circuit to execute all reporter-link changes together
        # Without recomputing the score change matrix.
        # Only actor-link changes require recomputation.
        score.change.matrix[a,b] = -Inf
        score.change.matrix[1:nActors,1:nActors] = -Inf
        max.change = max(score.change.matrix, na.rm = TRUE)
      }
    }
  }
}


#   Network score calculation
#   Assumes the following structure for lprobs:
#   lprobs$single.gt.wt is an n.actors x n.reporters matrix with the log probability
# that the single KO has a greater effect than the WT effect on the reporter.
#   lprobs$single.ngt.wt is log( 1 - exp( single.gt.wt ) )
#   lprobs$double.eq.single is a 3-d n.actors (1) x n.actors (2) x n.reporters array
# of the probability that actor 1 knocked out with actor 2 have the same effect on the
# reporter as actor 1 alone.
#   lprobs$double.gt.single is as above but probability of the double KO having a
# greater effect
scoreNetwork <- function( network, lp ){
  
  # Determine current ancestry.
  # ancestry(A,B) means A is an ancestor of B
  # start with edges
  ancestry <- network
  # grow the ancestry until no change
  repeat{
    old.ancestry <- ancestry
    for (node in 1:ncol(network)) {
      ancestors <- ancestry[,node]
      ancestry[,node] = ancestors | apply(network[,ancestors],1,any)
    }
    if( all(ancestry == old.ancestry) ) {break}
  }
  # Ok. now that we have the ancestry...
  
  # Score presence of an actor being an ancestor of a reporter
  reporter.ancestry <- ancestry[,howManyActors(lp)+(1:howManyReporters(lp))]
  score <-
    sum( ancestryScoreMatrix(lp)[reporter.ancestry], na.rm=TRUE ) +
    sum( log1mexp( ancestryScoreMatrix(lp)[!reporter.ancestry] ), na.rm=TRUE )
  
  #   Score two actors sharing a pathway: one is an ancestor of the other and both
  # are ancestors of the reporter.
  #   Score of two actors in independent pathways: both are ancestors of the reporter
  # but neither is an ancestor of the other.
  for (i1 in 2:howManyActors(lp)) {
    g1 = getActors(lp)[i1]
    for (i2 in 1:(i1-1)) {
      g2 = getActors(lp)[i2]
      
      both = reporter.ancestry[i1,] & reporter.ancestry[i2,]
      
      if( ancestry[i1,i2] || ancestry[i2,i1])
        score <- score + sum( scoreSharedPathways(lp, g1, g2)[both], na.rm=TRUE )
      else
        score <- score + sum( scoreIndependentPathways(lp, g1, g2)[both], na.rm=TRUE )
    }
  }
  
  return(score)
}

# Greedy reporter assignment:
# Finds the maximum-scoring network with arcs to reporters for a fixed actor network
getBestReporters = function(
  laps,
  nActors = howManyActors(laps),
  nReporters = howManyReporters(laps),
  network = matrix( FALSE, nrow=nActors, ncol = (nActors + nReporters) ) ){
  
  # We are not going to rely on scoreEdgeToggles for this (theoretically could?)
  # We are going to ignore the given score.
  
  # Debug statements!
  print("Getting best reporter configuration for the following actor network:")
  print(network[,1:nActors])
  
  # Get ancestry in actor network
  actorAncestry = inclusive.ancestry( network )
  rownames( actorAncestry ) =
    colnames( actorAncestry ) =
    rownames( network ) = getActors( laps )
  
  print("Computed ancestry:")
  print(actorAncestry)
  
  # Keep track of best contribution to score from each reporter
  scoreContributions = rep( -Inf, nReporters )
  
  # Keep track of seen ancestries
  ancestryHistory = matrix( nrow = nActors, ncol = 0)
  
  # For each subset of actors
  for ( n in 0:nActors ){
    for ( actors in combn( getActors( laps ), n, simplify = FALSE ) ){
      
      # Get full ancestry for this reporter combination
      ancestry = apply(as.matrix( actorAncestry[,actors] ),1,any)
      if (
        dim(ancestryHistory)[2] < 1 ||
        !any( apply( as.matrix( ancestry == ancestryHistory ), 2, all ) )
      ){
        ancestryHistory = cbind( ancestryHistory, ancestry )
        
        # For each reporter
        for ( r in 1:nReporters ){
          # Get the score contribution of this configuration:
          
          # Simple ancestry component:
          score =
            sum( ancestryScoreMatrix( laps )[ancestry,r], na.rm=TRUE ) +
            sum( log1mexp(ancestryScoreMatrix( laps )[!ancestry,r]), na.rm=TRUE )
          
          # 2le KO component:
          if ( length(actors) > 1 ){
            for ( i in 2:length(actors) ){
              for ( b in actors[1:(i-1)] ){
                a = actors[i]
                
                score = score +
                  if ( actorAncestry[a,b] || actorAncestry[b,a] )
                    scoreSharedPathways( laps, a, b )[r]
                  else
                    scoreIndependentPathways( laps, a, b )[r]
              }
            }
          }
          
          # Check against contribution vector
          if ( !is.na(score) && !is.nan(score) && score > scoreContributions[r] ){
            # Update contribution vector and network structure
            scoreContributions[r] = score
            network[,r] = FALSE
            network[actors,r] = TRUE
          }
        }
      }
    }
  }
  
  print( paste0( "Best score found to be: ", sum(scoreContributions) ) )
  
  infinite = which( scoreContributions == -Inf )
  if( length(infinite) > 0 ){
    print( paste0( "Found ", length(infinite), " prblematic reporters:" ) )
    print( getReporters( laps )[infinite] )
  }
  
  return( list( score = sum( scoreContributions ), network = network ) )
}

# Smart edge toggle
# laps stands for local ancestry and pathway scores
scoreEdgeToggles = function(
  network,
  laps) 
{
  # For easy writing
  n = howManyActors(laps)
  
  # The result matrix
  toggle.scores =
    matrix(
      nrow=n,
      ncol=n+howManyReporters(laps))
  
  actor.ancestry = inclusive.ancestry(network)
  
  # Get original ancestry vectors for reporters
  original.ancestry.vectors=matrix(nrow=n,ncol=howManyReporters(laps))
  for( r in 1:howManyReporters(laps))
    original.ancestry.vectors[,r] =
      apply(as.matrix(actor.ancestry[,network[,r+n]]),1,any)
  
  # Compute toggle scores for actor edges
  for( a in 1:n )
    for( b in (1:n)[-a] ){
      alternate.edges = network[1:n,1:n]
      alternate.edges[a,b] = !network[a,b]
      alternate.ancestry = inclusive.ancestry(alternate.edges)
      
      # score.change computation
      score.change = 0;
      for( r in 1:howManyReporters(laps) )
        score.change = score.change +
          compute.score.change(
            reporter = r,
            old.ancestry = original.ancestry.vectors[,r],
            new.ancestry = apply(as.matrix(alternate.ancestry[,network[,r+n]]),1,any),
            old.actor.ancestry = actor.ancestry,
            new.actor.ancestry = alternate.ancestry,
            laps,
            n)
      
      toggle.scores[a,b] = score.change
    }    
  
  # Compute toggle scores for reporter edges
  for( r in 1:howManyReporters(laps) )
    for( a in 1:n ){
      alternate.parent.vector = network[,r+n]
      alternate.parent.vector[a] = !alternate.parent.vector[a]
      alternate.ancestry.vector =
        apply(as.matrix(actor.ancestry[,alternate.parent.vector]),1,any)
      
      toggle.scores[a,r+n] = compute.score.change(
        reporter = r,
        old.ancestry = original.ancestry.vectors[,r],
        new.ancestry = alternate.ancestry.vector,
        old.actor.ancestry = actor.ancestry,
        new.actor.ancestry = actor.ancestry,
        laps,
        n)
    }

  return(toggle.scores)
}

compute.score.change = function(
  reporter,
  old.ancestry,
  new.ancestry,
  old.actor.ancestry,
  new.actor.ancestry,
  laps,
  n = howManyActors(laps))
{
  score.change = 0
  
  # Single-KO score update
  if( any( new.ancestry != old.ancestry ) )
    score.change = score.change +
      sum( ancestryScoreMatrix( laps )[new.ancestry,reporter], na.rm=TRUE ) +
      sum( log1mexp(ancestryScoreMatrix( laps )[!new.ancestry,reporter]), na.rm=TRUE ) -
      sum( ancestryScoreMatrix( laps )[old.ancestry,reporter], na.rm=TRUE ) -
      sum( log1mexp(ancestryScoreMatrix( laps )[!old.ancestry,reporter]), na.rm=TRUE )
  
  # Double-KO score update
  if( any(new.ancestry != old.ancestry) ||
      any(new.actor.ancestry != old.actor.ancestry, na.rm = TRUE))
    for( a in 2:n )
      for( b in 1:(a-1) ){
        
        if( new.ancestry[a] && new.ancestry[b] )
          if( new.actor.ancestry[a,b] || new.actor.ancestry[b,a] )
            score.change = score.change +
              scoreSharedPathways( laps, a, b )[reporter]
          else
            score.change = score.change +
              scoreIndependentPathways( laps, a, b )[reporter]
          
        if( old.ancestry[a] && old.ancestry[b] )
          if( old.actor.ancestry[a,b] || old.actor.ancestry[b,a] )
            score.change = score.change -
              scoreSharedPathways( laps, a, b )[reporter]
          else
            score.change = score.change -
              scoreIndependentPathways( laps, a, b )[reporter]
      }
  
  return(score.change)
}

# Inclusive ancestry of the square part of the matrix
inclusive.ancestry = function(network) {
  # Get square part of network
  rows = nrow( network )
  cols = ncol( network )
  n = min( rows, cols )
 
  # Add self-ancestry
  ancestry = network[1:n,1:n] | diag(n)
  repeat{
    old.ancestry = ancestry
    for (node in 1:n) {
      ancestry[,node] = apply(as.matrix(ancestry[,ancestry[,node]]),1,any)
    }
    if( all(ancestry == old.ancestry) ) { return(ancestry) }
  }
}

# Data prep help: get log probability from log-odds
lprob.from.lods = function( lods ){
  result = lods - log1p( exp( lods ) )
  result[ lods == Inf ] = 0
  return( result )
}

# And the backwards conversion
lprob2lods = function( lprob ){
  lprob - log1mexp( lprob )
}

# Data prep help: check that element-wise signs match on two vectors
signs.match = function( a, b ){ sign(a)==sign(b) }

# Compute log(1-exp(x)) accurately.
# Based on "Accurately Computing log(1-exp(-|a|))" by Martin Mächler
log1mexp = function(x){
  if ( length(x) < 1 ) return( numeric(0) )
  if ( length(x) > 1 ){
    result = log(-expm1(x))
    calculate.differently = is.finite(x) & (x < log(2))
    result[calculate.differently] = log1p(-exp(x[calculate.differently]))
    result
  } else {
    if( is.na(x) )
      NA
    else if( x < log(2) )
      log1p(-exp(x))
    else
      log(-expm1(x))
  }
}

# Utility for output
# Convert adjacency matrix to cytoscape edge list
edgesFromMatrix = function(
  adjacency,
  nodeNames =
    if(!is.null(colnames(adjacency)))
      colnames(adjacency)
    else
      1:max(dim(adjacency)) )
{
  indeces = which( adjacency, arr.ind=TRUE )
  
  sources = nodeNames[indeces[,1]]
  targets = nodeNames[indeces[,2]]

  # return
  data.frame( source=sources, target=targets, stringsAsFactors = FALSE )
}


# How to score an edge toggle
score.edge.toggles.old <- function( network, lprobs, score=score.network(network,lprobs) ){
  # Naive hackery:
  # Just score every move the old fashioned way
  toggle.scores <- matrix(0, nrow=nrow(network), ncol=ncol(network))
  for( i in 1:length(network) ){
    new.network <- network
    new.network[i] <- !network[i]
    toggle.scores[i] <- score.network(new.network,lprobs) - score
  }
  return(toggle.scores)
}