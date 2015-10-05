# A method for inferring the maximum-likelihood network based on log-odds from
# contrasts.
#
# Initial creation: August 10, 2015
# Yuriy Sverchkov

# Data structure holding relevant probabilities
# Might be used later.
#setClass("LogProbabilities", representation(
#  n.actors = "numeric", n.reporters = "numeric",
#  single.gt.wt = "numeric", double.gt.single = "numeric", double.eq.single = "numeric") )

# Constructs log-probability data structure or updates it.
# Semantics: for adding a single v. WT table use the KO gene for contrast.1
# For adding a double v. single table use the single KO gene for contrast.2 and the other for contrast.1
# Prior is the prior given to ebayes
logProbabilities <- function( lods, gt, prior, contrast.1, contrast.2 = "WT", logProbs = NULL ) {
  
  if( is.null(logProbs) ){
    logProbs$reporters = rownames( lods );
  }
  
  index.1 =
    if( contrast.1 %in% logProbs$actors )
      match(contrast.1, logProbs$actors)
    else{
      logProbs$actors = c( logProbs$actors, contrast.1 )
      length(logProbs$actors)
    }
  
  if( contrast.2 == "WT" ){# The single-KO tables
    
    lods = update.log.odds( lods, logProbs$single.gt.wt[index.1] - logProbs$single.ngt.wt[index.1], log(prior/(1-prior)))
    
    logProbs$single.gt.wt[index.1] = lprob.from.lods(lods) # TODO: add matching for differing reporter lists?
    logProbs$single.gt.wt[index.1][!gt] = -Inf
    logProbs$single.ngt.wt[index.1] = log1mexp( logProbs$single.ngt.wt )
  }else{# The double-KO tables
    index.2 =
      if( contrast.2 %in% logProbs$actors )
        match( contrast.2, logProbs$actors )
      else{
        logProbs$actors = c( logProbs$actors, contrast.2 )
        length(logProbs$actors)
      }
    
    lprobs = logProbs$double.eq.single[index.1][index.2]
    lods = update.log.odds( lods, log1mexp( lprobs ) - lprobs, log(prior/(1-prior)));
    
    lprobs = lprob.from.lods(lods)
    logProbs$double.eq.single[index.1][index.2] = log1mexp( lprobs )
    logProbs$double.gt.single[index.1][index.2] = lprobs
    logProbs$double.gt.single[index.1][index.2][~gt] = -Inf
  }

  logProbs
}

# Function for updating lods based on preexisting lods (log-prior-odds)
update.log.odds = function( log.odds.1, log.odds.2, log.prior.odds ){
  log.odds.1[is.na(log.odds.1)] = 0
  log.odds.2[is.na(log.odds.2)] = 0
  log.odds.1 + log.odds.2 - prior.odds
}

# Function for cleaning up log probs structure


# Multiple start (by default greedy) search:
multiStartNetworkSearch <- function(
  lprobs,
  n.actors=nrow(lprobs$single.gt.wt),
  n.reporters=ncol(lprobs$single.gt.wt),
  search.function=getMLNetwork)
  # See getMLNetwork for parameter definitions
{
  # We want to do a search starting with each combination of the possible
  # n.actors^2 - n.actors edges.
  # Start by getting the list of edges
  edges = which(
    matrix( TRUE, nrow=n.actors, ncol=n.actors ) &
    !diag( n.actors ) )
  
  # Keep track of max score
  max.score = -Inf
  max.network = NULL
  
  # Next iterate over 0-max num of edges
  for( n.edges in 0:length(edges) )
    # Iterate over combinations
    for( combination in combn(edges, n.edges, simplify = FALSE) ){
      # Make a starting network
      network = matrix( FALSE, nrow=n.actors, ncol=(n.actors+n.reporters) )
      for( i in combination ) network[i] = TRUE
      
      # Score
      candidate.score = score.network( network, lprobs )
      
      if( is.finite( candidate.score ) ){
        ml.result =
          search.function( lprobs, n.actors, n.reporters, network, candidate.score )
      
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

# Greedy edge toggle search:
getMLNetwork <- function(
  lprobs,
  n.actors=nrow(lprobs$single.gt.wt),
  n.reporters=ncol(lprobs$single.gt.wt),
  network=matrix( FALSE, nrow=n.actors,  ncol=(n.actors+n.reporters)),
  score=score.network( network, lprobs ) )
  # lods: data structure describing the log-odds comparing various KOs
  # n.actors: the number of things upstream in the network
  # n.reporters: the number of genes affected by the actors
  # network: initial network from which to begin the search (default is empty)
  # score: initial score of initial network
  # Graph is represented by the directed adjacency matrix where network[a,b] == TRUE
  # indicates the edge a->b
{
  repeat{
    score.change.matrix <- score.edge.toggles( network, lprobs )
    max.change <- max(score.change.matrix, na.rm = TRUE)
    if( max.change > 0 ){
      max.index <- which.max(score.change.matrix)
      network[max.index] = !network[max.index]
      score = score + max.change
    } else {
      return(list(network,score,score.change.matrix))
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
score.network <- function( network, lprobs, n.actors=nrow(network), n.reporters=(ncol(network)-nrow(network))){
  
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
  reporter.ancestry <- ancestry[,n.actors+(1:n.reporters)]
  score <-
    sum( lprobs$single.gt.wt[reporter.ancestry], na.rm=TRUE ) +
    sum( lprobs$single.ngt.wt[!reporter.ancestry], na.rm=TRUE )
  
  #   Score two actors sharing a pathway: one is an ancestor of the other and both
  # are ancestors of the reporter.
  #   Score of two actors in independent pathways: both are ancestors of the reporter
  # but neither is an ancestor of the other.
  for (g1 in 2:n.actors) {
    for (g2 in 1:(g1-1)) {
      both = reporter.ancestry[g1,] & reporter.ancestry[g2,]
      if( ancestry[g1,g2] || ancestry[g2,g1]){
        score <- score +
          sum( lprobs$double.eq.single[g1,g2,both], na.rm=TRUE ) +
          sum( lprobs$double.eq.single[g2,g1,both], na.rm=TRUE )# +
          #sum( logspace.not(lprobs$double.eq.single[g1,g2,!both]), na.rm=TRUE ) +
          #sum( logspace.not(lprobs$double.eq.single[g2,g1,!both]), na.rm=TRUE )
      }else{
        score <- score +
          sum( lprobs$double.gt.single[g1,g2,both], na.rm=TRUE ) +
          sum( lprobs$double.gt.single[g2,g1,both], na.rm=TRUE )# +
          #sum( logspace.not(lprobs$double.gt.single[g1,g2,!both]), na.rm=TRUE ) +
          #sum( logspace.not(lprobs$double.gt.single[g2,g1,!both]), na.rm=TRUE )
      }
    }
  }
  
  return(score)
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

# Smart edge toggle
score.edge.toggles = function(
  network,
  local.scores,
  actor.ancestry = inclusive.ancestry(network),
  n.actors = nrow(network),
  n.reporters = ncol(network) - n.actors)
{
  # The result matrix
  toggle.scores = matrix(nrow=n.actors, ncol=n.actors+n.reporters)
  
  # Get original ancestry vectors for reporters
  original.ancestry.vectors=matrix(nrow=n.actors,ncol=n.reporters)
  for( r in 1:n.reporters)
    original.ancestry.vectors[,r] =
      apply(as.matrix(actor.ancestry[,network[,r+n.actors]]),1,any)
  
  # Compute toggle scores for actor edges
  for( a in 1:n.actors )
    for( b in (1:n.actors)[-a] ){
      alternate.edges = network[1:n.actors,1:n.actors]
      alternate.edges[a,b] = !network[a,b]
      alternate.ancestry = inclusive.ancestry(alternate.edges)
      
      # score.change computation
      score.change = 0;
      for( r in 1:n.reporters )
        score.change = score.change +
          compute.score.change(
            reporter = r,
            old.ancestry = original.ancestry.vectors[,r],
            new.ancestry = apply(as.matrix(alternate.ancestry[,network[,r+n.actors]]),1,any),
            old.actor.ancestry = actor.ancestry,
            new.actor.ancestry = alternate.ancestry,
            local.scores,
            n.actors)
      
      toggle.scores[a,b] = score.change
    }    
  
  # Compute toggle scores for reporter edges
  for( r in 1:n.reporters )
    for( a in 1:n.actors ){
      alternate.parent.vector = network[,r+n.actors]
      alternate.parent.vector[a] = !alternate.parent.vector[a]
      alternate.ancestry.vector =
        apply(as.matrix(actor.ancestry[,alternate.parent.vector]),1,any)
      
      toggle.scores[a,r+n.actors] = compute.score.change(
        reporter = r,
        old.ancestry = original.ancestry.vectors[,r],
        new.ancestry = alternate.ancestry.vector,
        old.actor.ancestry = actor.ancestry,
        new.actor.ancestry = actor.ancestry,
        local.scores,
        n.actors)
    }

  return(toggle.scores)
}

compute.score.change = function(
  reporter,
  old.ancestry,
  new.ancestry,
  old.actor.ancestry,
  new.actor.ancestry,
  local.scores,
  n.actors = nrows(actor.ancestry))
{
  score.change = 0
  
  # Single-KO score update
  if( any( new.ancestry != old.ancestry ) )
    score.change = score.change +
      sum( local.scores$single.gt.wt[new.ancestry,reporter], na.rm=TRUE ) +
      sum( local.scores$single.ngt.wt[!new.ancestry,reporter], na.rm=TRUE ) -
      sum( local.scores$single.gt.wt[old.ancestry,reporter], na.rm=TRUE ) -
      sum( local.scores$single.ngt.wt[!old.ancestry,reporter], na.rm=TRUE )
  
  # Double-KO score update
  if( any(new.ancestry != old.ancestry) ||
      any(new.actor.ancestry != old.actor.ancestry, na.rm = TRUE))
    for( a in 2:n.actors )
      for( b in 1:(a-1) ){
        
        if( new.ancestry[a] && new.ancestry[b] )
          if( new.actor.ancestry[a,b] || new.actor.ancestry[b,a] )
            score.change = score.change +
              sum( local.scores$double.eq.single[a,b,reporter], na.rm=TRUE ) +
              sum( local.scores$double.eq.single[b,a,reporter], na.rm=TRUE )
          else
            score.change = score.change +
              sum( local.scores$double.gt.single[a,b,reporter], na.rm=TRUE ) +
              sum( local.scores$double.gt.single[b,a,reporter], na.rm=TRUE )
          
        if( old.ancestry[a] && old.ancestry[b] )
          if( old.actor.ancestry[a,b] || old.actor.ancestry[b,a] )
            score.change = score.change -
              sum( local.scores$double.eq.single[a,b,reporter], na.rm=TRUE ) -
              sum( local.scores$double.eq.single[b,a,reporter], na.rm=TRUE )
          else
            score.change = score.change -
              sum( local.scores$double.gt.single[a,b,reporter], na.rm=TRUE ) -
              sum( local.scores$double.gt.single[b,a,reporter], na.rm=TRUE )
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
  return( lods - log1p(1+expm1( lods )) )
}

# Data prep help: check that element-wise signs match on two vectors
signs.match = function( a, b ){ sign(a)==sign(b) }

# Compute log(1-exp(x)) accurately.
# Based on "Accurately Computing log(1-exp(-|a|))" by Martin Mächler
log1mexp = function(x){
  if( length(x) > 1 ){
    result = log(-expm1(x))
    calculate.differently = is.finite(x) & (x < log(2))
    result[calculate.differently] = log1p(-exp(x[calculate.differently]))
    result
  }else{
    if( x < log(2) )
      log1p(-exp(x))
    else
      log(-expm1(x))
  }
}

# Utility for output
# Convert adjacency matrix to cytoscape edge list
edges.from.matrix = function(
  adjacency.matrix,
  node.names =
    if(!is.null(colnames(adjacency.matrix)))
      colnames(adjacency.matrix)
    else
      1:max(dim(adjacency.matrix)) )
{
  n = sum(adjacency.matrix);
  edges = data.frame( source=rep("",n), target=rep("",n), stringsAsFactors = FALSE )
  edge = 1
  for( i in 1:nrow(adjacency.matrix) )
    for( j in 1:ncol(adjacency.matrix) )
      if( adjacency.matrix[i,j]){
        edges[edge,"source"] = node.names[i];
        edges[edge,"target"] = node.names[j];
        edge = edge+1;
      }
  
  return(edges);
}