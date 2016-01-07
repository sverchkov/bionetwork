## Functions having to do with ancestry aggregation etc.
##

#' Aggregate individual reporter-ancestries into one master ancestry
#' 
#' @param ancestries - a actors x actors+1 x reporters array of individual
#' reporter-ancestries.
#' @return a single ancestry (with possibly duplicated actors) which is a matrix of
#' size (new actors) x (new actors + reporters)
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
  
  # An indexing system to keep track of actor copies
  actor.versions = list();
  for ( i in 1:nA ){
    actor.versions[[ actors[i] ]]$matrix = matrix( nrow = nA, ncol = 0, dimnames = list( actors, vector() ) )
    actor.versions[[ actors[i] ]]$indeces = vector( mode = "numeric" )
  }

  # First pass to establish the actor set
  numberOfActorNodes = 0
  
  for ( reporter in reporters ) { # For each reporter
    # Get its ancestors
    ancestry = ancestries[, nA+1, reporter ]
    
    for ( actor in actors( ancestry ) ) { # For each ancestor
      
      actor.ancestors = ancestries[, actor, reporter]
      
      # Check if it is already in the version list
      version.matrix = actor.versions[[ actor ]]$matrix
      v = which( apply( version.matrix == actor.ancestors, 2, all ) )
      if ( v == 0 ) { # If not, add it
        v = dim( version.matrix )[ 2 ] + 1
        version.matrix[, v ] = actor.ancestors
        actor.versions[[ actor ]]$matrix = version.matrix
        numberOfActorNodes = numberOfActorNodes + 1
        actor.versions[[ actor ]]$indeces[v] = numberOfActorNodes
      }
    }
  }
  
  # Actor node name vector
  actor.node.names = vector( mode = "character", length = numberOfActorNodes )
  for ( actor in actors ) {
    indeces = actor.versions[[ actor ]]$indeces;
    for ( i in 1:length( indeces ) ) {
      actor.node.names[ indeces[ i ] ] = paste0( actor, ".", i )
    }
  }
  # Result ancestry, now with new actor nodes.
  result.ancestry = matrix(
    data = FALSE,
    nrow = numberOfActorNodes,
    ncol = numberOfActorNodes + nR,
    dimnames = list(
      actor.node.names,
      c( actor.node.names, reporters )
    ) )
    
  # Second pass builds the result matrix
  for ( reporter in reporters ) { # For each reporter
    # Get its ancestors
    ancestry = ancestries[, nA+1, reporter ]
    
    for ( actor in actors( ancestry ) ) { # For each ancestor
      
      actor.ancestors = ancestries[, actor, reporter]
      
      # Get the actor index
      i = actor.versions[[ actor ]]$indeces[ which( apply( actor.versions[[ actor ]]$matrix == actor.ancestors, 2, all ) ) ]
      
      result.ancestry[ i, reporter ] = TRUE
      
      for ( ancestor in actors( actor.ancestors ) ){
        
        # Get the ancestor index
        j = actor.versions[[ ancestor ]]$indeces[ which( apply( actor.versions[[ ancestor ]]$matrix == ancestries[, ancestor, reporter], 2, all ) ) ]
        result.ancestry[ j, i ] = TRUE
      }
    }
  }
  
  result.ancestry
}