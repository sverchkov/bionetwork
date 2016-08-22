# Graph scorer, version concieved August 18, 2016
# By Yuriy Sverchkov

#' Get the effect-speciic score for graph adj and effect link vector theta for a specific effect
#' 
#' @param adj The AxA adjacency matrix. Should have 1s on the diagonal.
#' @param theta The length-A effect link vector, 1 for direct links to effect from action.
#' @param vectorSA The length-A vector of single-action log-odds that action induced a change in the effect
#' @param matrixDA The AxA matrix that the double-action i,j is differentiall expressed from action i alone
#' @return The score
#' @export
score4effect = function ( adj, theta, vectorSA, matrixDA ){
  # Predicted affecting actions vector:
  predictedSA = adj %*% theta > 0
  
  # Predicted affecting
  term2 = !adj & !sweep( adj, 2, theta, `&` )
  predictedDA = sweep( term2, 1, predictedSA, `&` )
  
  return ( sum( predictedSA * vectorSA ) + sum( predictedDA * matrixDA ) )
}

#' Get the best theta and effect-specific score for graph adj.
#' 
#' @param adj The AxA adacency matrix. Should have 1s on the diagonal.
#' @param vectorSA The length-A vector of single-action log-odds that action induced a change in the effect
#' @param matrixDA The AxA matrix that the double-action i,j is differentiall expressed from action i alone
#' @param regularization A regularization constant that penalizes nonzero theta-elements. Default is 0 (no regularization)
#' @return The best theta-assignment and score, in a named list
#' @export
scoreBestTheta = function ( adj, vectorSA, matrixDA, regularization = 0 ){
  best.score = -Inf
  best.theta = NULL
  n = dim( adj )[1]
  if ( n >= 31 ) stop( "Sorry, this isn't written to handle more than 30 actions :(" )
  
  for ( i in 1:(2^n) ) {
    theta = as.logical( intToBits( i )[1:n] )
    #print( theta )
    score = ( 1 - regularization ) * score4effect( adj, theta, vectorSA, matrixDA ) - regularization * sum( theta )
    #print( c( score, best.score ) )
    if ( score > best.score ) {
      best.score = score
      best.theta = theta
    }
  }
  
  return ( list( score = best.score, theta = best.theta ) )
}

#' Get the MAP of thetas for a specific graph
#' 
#' @param adj the AxA adjacency matrix. Should have 1s on the diagonal
#' @param matrixSA The AxE matrix of log odds that action a induced a change in the effect e
#' @param tensorDA An AxAxE array where element [a,b,e] is the log-odds that effect e in DKO ab is DE W.R.T. SKO a
#' @param regularization A regularization constant between 0 and 1 that penalizes nonzero theta-elements
#' @return The best theta-assignment for each effect and the associated scores, in a named list
#' @export
scoreMAPTheta = function ( adj, matrixSA, tensorDA, regularization = 0 ){
  
  actions = rownames( matrixSA )
  effects = colnames( matrixSA )
  
  n.actions = length( actions )
  n.effects = length( effects )
  
  best.scores = rep( -Inf, n.effects )
  best.theta = matrix( FALSE, n.actions, n.effects, dimnames = list( actions, effects ) )
  
  tensor = tensorDA
  for ( action in actions )
    tensor[ action, action, ] = matrixSA[ action, ]
  
  if ( n.actions >= 31 ) stop( "Sorry, this isn't written to hangle more than 30 actions :(" )
  
  for ( i in 1:(2^n.actions) ) {
    theta = as.logical( intToBits( i )[ 1:n.actions ] )
    scores =
      ( 1 - regularization ) * apply( tensor, 3, function ( m ) score4effect( adj, theta, diag( m ), m ) ) -
      regularization * sum( theta )
    
    better = ( scores > best.scores )
    best.theta[ , better ] = theta
    best.scores[ better ] = scores[ better ]
  }
  
  return ( list( score = best.scores, theta = best.theta ) )
}