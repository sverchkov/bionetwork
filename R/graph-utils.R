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
      ancestor.row.select = ancestry[,j]
      if ( n < m ) ancestor.col.select = c( ancestor.col.select, rep( FALSE, m-n ) )
      if ( n > m ) ancestor.col.select = ancestor.col.select[1:m]
      ancestry[,j] = apply( ancestry[,ancestor.col.select], 1, any ) & ancestry[,j]
    }
  }
  
  # return
  ancestry
}