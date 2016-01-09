## Functions having to do with ancestry aggregation etc.
##

#' Aggregate individual reporter-ancestries into one master ancestry
#' 
#' @param ancestries - a actors x actors+1 x reporters array of individual
#' reporter-ancestries.
#' @return a single ancestry (with possibly duplicated actors) which is a matrix of
#' size (new actors) x (new actors + reporters)
#' @export
aggregateAncestries = function ( ancestries ){
  
  # Get dimensions and names
  ds = dim( ancestries )
  nA = ds[1]
  nR = ds[3]
  
  dn = dimnames( ancestries )
  
  actorMap = vector( mode = "character", length = 0 )
  actors = 1:nA
  reporters = 1:nR
  
  if ( ! is.null( dn ) ){
    actors = dn[[ 1 ]]
    reporters = dn[[ 3 ]]
  }  
  
  # An indexing system to keep track of actor copies
  actor.versions = list();
  for ( i in 1:nA ){
    actor.versions[[ actors[i] ]]$array = array( dim = c( nA, nA, 0 ), dimnames = list( actors, actors, vector() ) )
    actor.versions[[ actors[i] ]]$indeces = vector( mode = "numeric" )
  }

  # First pass to establish the actor set
  numberOfActorNodes = 0
  
  for ( reporter in reporters ) { # For each reporter
    # Get its ancestors
    ancestry = ancestries[, nA+1, reporter ]
    
    for ( actor in actors[ ancestry ] ) { # For each ancestor
      
      actor.ancestors = ancestries[, 1:nA, reporter]
      actor.ancestors[, ! ancestries[, actor, reporter ] ] = FALSE
      
      # Add if it is already in the version list
      version.matrix = actor.versions[[ actor ]]$array

      # The c around actor.anestors is a workaround around a non-comformable-arrays error      
      if ( length( version.matrix ) < 1 || !any( apply( version.matrix == c( actor.ancestors ), 3, all ) ) ) {
        version.matrix = abind::abind( version.matrix, actor.ancestors, along = 3 )
        actor.versions[[ actor ]]$array = version.matrix
        numberOfActorNodes = numberOfActorNodes + 1
        actor.versions[[ actor ]]$indeces = c( actor.versions[[ actor ]]$indeces, numberOfActorNodes )
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
    
    for ( actor in actors[ ancestry ] ) { # For each ancestor
      
      actor.ancestors = ancestries[, 1:nA, reporter]
      actor.ancestors[, ! ancestries[, actor, reporter ] ] = FALSE
      
      # Get the actor index
      i = actor.versions[[ actor ]]$indeces[ which( apply( actor.versions[[ actor ]]$array == c(actor.ancestors), 3, all ) ) ]
      
      result.ancestry[ i, reporter ] = TRUE
      
<<<<<<< HEAD
      for ( ancestor in actors[ actor.ancestors ] ){
=======
      for ( ancestor in actors[ ancestries[, actor, reporter ] ] ){

        ancestor.ancestors = ancestries[, 1:nA, reporter]
        ancestor.ancestors[, ! ancestries[, ancestor, reporter ] ] = FALSE
>>>>>>> 3edfe453cc2332941f1e250dd48158b9097d8ea7
        
        # Get the ancestor index
        j = actor.versions[[ ancestor ]]$indeces[ which( apply( actor.versions[[ ancestor ]]$array == c(ancestor.ancestors), 3, all ) ) ]
        result.ancestry[ j, i ] = TRUE
      }
    }
  }
  
  result.ancestry
}