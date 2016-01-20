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
  d = min( n, m ) # Diagonal length
  
  ancestry = adjacency # Parents are ancestors
  old_ancestry = matrix( FALSE, n, m )
  
  while ( any( ancestry != old_ancestry ) ){
    old_ancestry = ancestry
    for ( j in 1:m ){
      anc = as.matrix( ancestry[ which( ancestry[,j] ) ] )
      if ( any( anc ) ) {
        ancestry[,j] = apply( anc, 1, any ) & ancestry[,j]
      }
    }
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