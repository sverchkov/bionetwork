# Miscellaneaous utilities

#' Make a contrast matrix
mkContrastMatrix = function( single.genes, double.specs, wt.str, the.colnames ){
  
  n.1le = length(single.genes)
  n.2le = length(double.specs)
  n.col = n.1le+2*n.2le
  
  # Contrast generation:
  # We make one contrast for each singke KO and two contrasts for each double KO.
  contrast.matrix = matrix( data = 0, nrow = length(the.colnames), ncol = n.col )
  
  col.names = c( single.genes, vector(mode = "character", length = 2*n.2le ) )
  
  # Single KO contrasts:
  i = 0
  for( gene in single.genes ){
    i = i+1
    contrast.matrix[which(the.colnames == gene),i] = 1
    contrast.matrix[which(the.colnames == wt.str),i] = -1
  }
  # Double KO contrasts:
  for( element in double.specs )
    for( single in element[[2]] ){
      i = i+1
      contrast.matrix[which(the.colnames == element[[1]]),i] = 1
      contrast.matrix[which(the.colnames == single),i] = -1
      col.names[i] = paste0(element[[1]],'-',single)
    }
  
  colnames( contrast.matrix ) = col.names
  rownames( contrast.matrix ) = the.colnames
  
  contrast.matrix
}
