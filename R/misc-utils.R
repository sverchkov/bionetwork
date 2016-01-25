# Miscellaneaous utilities

#' Make a contrast matrix
mkContrastMatrix = function( single.genes, double.specs, wt.str, the.colnames ){
  
  n.1le = length(single.genes)
  n.2le = length(double.specs)
  n.col = n.1le+3*n.2le
  
  # Contrast generation:
  # We make one contrast for each singke KO and two contrasts for each double KO.
  contrast.matrix = matrix( data = 0, nrow = length(the.colnames), ncol = n.col )
  
  col.names = c( single.genes, vector(mode = "character", length = 3*n.2le ) )
  
  # Single KO contrasts:
  i = 0
  for( gene in single.genes ){
    i = i+1
    contrast.matrix[which(the.colnames == gene),i] = 1
    contrast.matrix[which(the.colnames == wt.str),i] = -1
  }
  # Double KO contrasts:
  for( element in double.specs )
    for( single in c( element[[2]], wt.str ) ){
      i = i+1
      contrast.matrix[which(the.colnames == element[[1]]),i] = 1
      contrast.matrix[which(the.colnames == single),i] = -1
      col.names[i] = paste0(element[[1]],'-',single)
    }
  
  colnames( contrast.matrix ) = col.names
  rownames( contrast.matrix ) = the.colnames
  
  contrast.matrix
}

#' Get groups of reporters based on parent signatures.
groupReporters = function( network, name.groups = TRUE ){
  # network is taken to be the (actors x actors+reporters) adjacency matrix, with column
  # names annotated.
  # name.groups = TRUE uses the network's row names to construct a space-separated list
  # for each parent set. When to FALSE, it returns numbered groups, the index-1 being the
  # bitmap defining the parent set.
  
  d = dim( network )
  n = d[1]     # n = number of actors
  m = d[2] - n # m = number of reporters
  
  nGroups = 2^n
  groups = list()
  for ( i in 1:nGroups ){
    parents = ( 1 == intToBits( i-1 )[1:n] )
    group = names( which( apply( parents == network[, n + (1:m) ], 2, all ) ) )
    if ( name.groups ){
      if ( length( group ) > 0 ){
        groupName = paste0("{",paste(rownames(network)[parents], collapse=" "),"}")
        groups[[groupName]] = group
      }
    } else {
      groups[[i]] = group
    }
  }
  return( groups )
}

#' Print groups to csv
list2csv = function( lst, file ){
  sink( file )
  for ( name in names( lst ) ){
    cat( name )
    cat( ", " )
    cat( paste( lst[[name]], collapse = ", " ) )
    cat( "\n" )
  }
  sink()
}

#' Replaces all NAs in any indexable data structure (with 0 by default)
#' 
#' @param x - original
#' @param replacement - the thing with which to replace zeros
#' @return version where all NAs were replaced by 0s
replaceNAs = function ( x, replacement = 0 ) {
  x[ is.na( x ) ] = replacement
  return ( x )
}

#' Make a logic vector
#' 
#' Makes a vector of a cetain length with a true value in a certain position (remainder false)
#' @param n vector length
#' @param i the position of the true value
logicVector = function( n = 1, i = 1 ) c(rep(FALSE,i-1),TRUE,rep(FALSE,n-i))