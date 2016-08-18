# Graph utility functions

#' Converts an adjacency matrix to an ancestry matrix
#'
#' Converts an adjacency matrix (where element [i,j] is true when there is an edge from i to j) to an ancestry matrix
#' (where element [i,j] is true when i is an ancestor of j).
#' This function follows the convention that every node is considered its own ancestor.
#' @param adjacency - the n x m adjacency matrix
#' @return an n x m ancestry matrix
adjacencyToAncestry = function ( adjacency ){
  n = nrow( adjacency ) # Number of rows
  m = ncol( adjacency ) # Number of columns
  if ( n > m )
    stop( "Can't convert adjacency with rows > cols to ancestry!")
  
  ancestry = adjacency # Parents are ancestors
  old_ancestry = matrix( FALSE, n, m )
  
  while ( any( ancestry != old_ancestry ) ){
    old_ancestry = ancestry

    ancestry = ancestry | ( ancestry[1:n,1:n] %*% ancestry )
    # OLD WAY
    #for ( j in 1:m ){
    #  anc = as.matrix( ancestry[ , which( ancesmatritry[,j] ) ] )
    #  if ( any( anc ) ) {
    #    ancestry[,j] = apply( anc, 1, any ) & ancestry[,j]
    #  }
    #}
  }
  
  # return
  ancestry
}

#' Converts an ancestry matrix to an adjacency matrix
#' 
#' Converts an ancestry matrix (where element [i,j] is true when i is an ancestor of j)
#' to an adjacency matrix (where element [i,j] is true when there is an edge from i to j).
#' This function does not take self-ancestry seriously.
#' @param ancestry - the n x m ancestry matrix
#' @return an n x m adjacency matrix
ancestryToAdjacency = function ( ancestry ){
  n = nrow( ancestry )
  m = ncol( ancestry )
  d = min( n, m ) # Diagonal length
  
  stop("Not yet implemented")
  
  adjacency = ancestry & (0 == diag( nrow = n, ncol = m ) )
  
  # Go through nodes, and clear their ancestors
  for ( i in 1:m ){
    
    descendants = adjacency[i,]
    adjacency[, descendants ] = adjacency[, descendants ] & !adjacency[, i ]
  }
  
  return ( adjacency )
}

#' Transitively closed adjacency matrix
#' 
#' Returns true if adjacency matrix is transitively closed, note that self-edges count.
#' @param adj The adjacency matrix
#' @return True iff the adjacency matrix is transitively closed
transitivelyClosed = function ( adj ){
  all( adj == ( ( adj %*% adj ) > 0 ) )
}

#' Enumerate all matrices that have more edges
#'
#' @param adj The adjacency matrix
#' @param firstEdge First edge to consider adding ( default = 1 )
#' @return The list of matrices obtained from adding edges to this one ( including input matrix )
allFullerMatrices = function ( adj, firstEdge = 1 ){

  edges = which( !adj )
  edge = min( edges[ edges >= firstEdge ] )
  if ( is.finite( edge ) ){
    adj2 = adj
    adj2[ edge ] = TRUE
    return ( c( allFullerMatrices( adj, edge+1 ), allFullerMatrices( adj2, edge+1 ) ) )
  } else return ( list( adj ) )
}

#' Enumerate all transitively closed n x n matrices
#' 
#' @param n The number of rows and columns
#' @param fullDiagonal Whether to enforce 1's on the diagonal (defult = TRUE)
#' @return All transitively closed n x n matrices
enumerateTransitivelyClosedMatrices = function ( n, fullDiagonal = TRUE ) {
  M = if ( fullDiagonal ) diag( n ) > 0 else matrix( TRUE, n, n )
  matrices = allFullerMatrices( M )
  return ( matrices[ sapply( matrices, transitivelyClosed ) ] )
}