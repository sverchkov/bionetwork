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

# Print groups to csv
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

## OLD STUFF. NOT SURE IF ANY OF THIS IS USED:
# 
# special.order <- function( matrix, column.list ){
#   if( length( column.list ) > 1 ){
#     matrix <- special.order( matrix, column.list[-1] )
#   }
#   if( length( column.list ) > 0 ){
#     matrix <- matrix[order(abs(matrix[,column.list[1]]),matrix[,column.list[1]],decreasing = TRUE),]
#   }
#   return(matrix)
# }
# 
# select.suspicious.genes <- function( matrix, reference.column, columns ){
#   matches = FALSE
#   for( c in columns ){
#     matches = matches | (matrix[,reference.column]*matrix[,c]==-1)
#   }
#   return(matrix[matches,])
# }
# 
# select.eq.to.ref <- function( matrix, reference.column, columns ){
#   matches = FALSE
#   for( c in columns ){
#     matches = matches | (matrix[,reference.column] == matrix[,c])
#   }
#   return(matrix[matches,])
# }