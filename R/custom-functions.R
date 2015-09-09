special.order <- function( matrix, column.list ){
  if( length( column.list ) > 1 ){
    matrix <- special.order( matrix, column.list[-1] )
  }
  if( length( column.list ) > 0 ){
    matrix <- matrix[order(abs(matrix[,column.list[1]]),matrix[,column.list[1]],decreasing = TRUE),]
  }
  return(matrix)
}

select.suspicious.genes <- function( matrix, reference.column, columns ){
  matches = FALSE
  for( c in columns ){
    matches = matches | (matrix[,reference.column]*matrix[,c]==-1)
  }
  return(matrix[matches,])
}

select.eq.to.ref <- function( matrix, reference.column, columns ){
  matches = FALSE
  for( c in columns ){
    matches = matches | (matrix[,reference.column] == matrix[,c])
  }
  return(matrix[matches,])
}