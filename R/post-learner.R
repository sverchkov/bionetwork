# Relearning with penalty
# To resolve never-ending computations, we take the multi-pool network structure and
# re-attach all the reporters to it.
#
# To avoid overfitting, we also repeat the process with certain duplicate actors
# removed.

#' Greedy reporter assignment
#' 
#' Finds the maximum-scoring network with arcs to reporters for a fixed actor network
#' using ancestries.
#' @param lll The LocalLogLikelihoods object
#' @param actor.ancestry The actor ancestry
#' @param debug.statements Whetehr debug statements should be printed
#' @return an ancestry where all reporters are attached and the associated score
#' @export
attachReporters = function(
  lll,
  actor.ancestry,
  debug.statements = FALSE ){
  
  # Number of actor versions
  nA = nrow( actor.ancestry )
  
  # Number of reporters
  nR = howManyReporters( lll )
  
  # Actor list
  actor.list = rownames( actor.ancestry )
  
  # Ensure square ancestry.
  actor.ancestry = actor.ancestry[1:nA,1:nA]
  
  if ( debug.statements ){# Debug statements!
    print("Getting best reporter configuration for the following actor network:")
    print( actor.ancestry )
  }
  
  # Keep track of best contribution to score from each reporter
  score.contributions = rep( -Inf, nR )
  reporter.ancestries = matrix( FALSE, nrow = nA, ncol = nR, dimnames = list( actor.list, getReporters( lll ) ) )
  
  # Keep track of seen ancestries
  ancestry.history = matrix( nrow = nA, ncol = 0)
  
  # For each subset of actors
  for ( n in 0:nA ){
    for ( actors in combn( actor.list, n, simplify = FALSE ) ){
      
      # Get full ancestry for this reporter combination
      ancestry = apply( as.matrix( actor.ancestry[, actors ] ), 1, any )
      if (
        !any( duplicated( stripDots( actor.list[ ancestry ] ) ) ) &&
        ( dim(ancestry.history)[2] < 1 ||
        !any( apply( as.matrix( ancestry == ancestry.history ), 2, all ) ) )
      ){
        ancestry.history = cbind( ancestry.history, ancestry )
        
        # Get log-likelihood contribution for each reporter
        contributions = scoreLikelihoodsPerReporter(
          cbind( actor.ancestry,
                 matrix( rep( ancestry, nR ), nA, nR,
                         dimnames =
                           list( actor.list,
                                 getReporters( lll ) )
                         )
                 ), lll )
        
        winners = contributions > score.contributions
        if( any( winners ) ){
          score.contributions[winners] = contributions[winners]
          reporter.ancestries[,winners] = ancestry
        }
      }
    }
  }
  
  if ( debug.statements ){
    print( paste0( "Best score found to be: ", sum( score.contributions ) ) )
  
    infinite = which( score.contributions == -Inf )
    if( length( infinite ) > 0 ){
      print( paste0( "Found ", length(infinite), " prblematic reporters:" ) )
      print( getReporters( laps )[infinite] )
    }
  }
  
  return( list( score = sum( score.contributions ), ancestry = cbind( actor.ancestry, reporter.ancestries ) ) )
}

#' Strip dots from strings
#' 
#' Only keeps the part of each string that preceeds the first dot (.) character.
#' @param strings - the original string array
#' @return strings with everything after the first dot (inclusive) stripped
stripDots = function ( strings ){
  sapply( strsplit( strings, ".", fixed = TRUE ), function ( x ) x[1] )
}