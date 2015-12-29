# This is the new likelihood-based scoring class and its functions
#' @import limma
#' @export LocalLogLikelihoods makeLLL

# Class definition
LocalLogLikelihoods = setClass( "LocalLogLikelihoods",
                                slots = c( actors = "character"
                                         , reporters = "character"
                                         , eqLogLik = "matrix" # Log-likelihoods given equal expression change
                                         , difLogLik = "matrix" # Log-likelihoods given different expression change
                                         , wt = "character" # The wild-type string
                                        #, ampLogLik = "matrix" # Log-likelihoods given amplified expression change
                                ) )

# Simple property functions
setMethod(f = "getActors",
          signature = "LocalLogLikelihoods",
          definition = function(theObject) theObject@actors )

setMethod(f = "getReporters",
          signature = "LocalLogLikelihoods",
          definition = function ( theObject ) theObject@reporters )

# Interfaces for scoring function
setMethod(f = "nonAncestryScoreMatrix",
          signature = "LocalLogLikelihoods",
          definition = function ( theObject ) t( theObject@eqLogLik[, getActors( theObject ) ] ) )

setMethod(f = "ancestryScoreMatrix",
          signature = "LocalLogLikelihoods",
          definition = function ( theObject ) t( theObject@difLogLik[, getActors( theObject ) ] ) )

setMethod(f = "scoreSharedPathways",
          signature( theObject = "LocalLogLikelihoods", actor1 = "character", actor2 = "character" ),
          definition = function ( theObject, actor1, actor2 ){
            pmax(
              accessTable( theObject@eqLogLik, actor1, actor2, actor1 ),
              accessTable( theObject@eqLogLik, actor1, actor2, actor2 ) )
          } )

setMethod(f = "scoreIndependentPathways",
          signature( theObject = "LocalLogLikelihoods", actor1 = "character", actor2 = "character" ),
          definition = function ( theObject, actor1, actor2 ){
            pmin(
              accessTable( theObject@difLogLik, actor1, actor2, actor1 ),
              accessTable( theObject@difLogLik, actor1, actor2, actor2 )
            )
          } )

setMethod( f = "scoreSingleAncestor",
           signature( theObject = "LocalLogLikelihoods", ancestor = "character", nonancestor = "character" ),
           definition = function ( theObject, ancestor, nonancestor ) {
             accessTable( theObject@eqLogLik, ancestor, nonancestor, ancestor )
           } )

setMethod( f = "scoreNeitherAncestor",
           signature( theObject = "LocalLogLikelihoods", actor1 = "character", actor2 = "character" ),
           definition = function ( theObject, actor1, actor2 ) {
             accessTable( theObject@eqLogLik, actor1, actor2, theObject@wt )
           } )

# We might write in numeric equivalents later. For now the char enforcement should do.

#' Overloaded + operator to merge LocalLogLikelihoods objects
setMethod( "+", signature( e1 = "LocalLogLikelihoods", e2 = "LocalLogLikelihoods" ), function ( e1, e2 ){
  
  reporters = union( e1@reporters, e2@reporters )
  actors = union( e1@actors, e2@actors )
  columns1 = colnames( e1@eqLogLik )
  columns2 = colnames( e2@eqLogLik )
  table.columns = union( columns1, columns2 )
  
  difLogLik = eqLogLik = matrix( data = 0,
                                 nrow = length( reporters ),
                                 ncol = length( table.columns ),
                                 dimnames = list( reporters, table.columns ))
  
  difLogLik[ e1@reporters, columns1 ] = difLogLik[ e1@reporters, columns1 ] + e1@difLogLik[ e1@reporters, columns1 ]
  difLogLik[ e2@reporters, columns2 ] = difLogLik[ e2@reporters, columns2 ] + e2@difLogLik[ e2@reporters, columns2 ]
  eqLogLik[ e1@reporters, columns1 ] = eqLogLik[ e1@reporters, columns1 ] + e1@eqLogLik[ e1@reporters, columns1 ]
  eqLogLik[ e2@reporters, columns2 ] = eqLogLik[ e2@reporters, columns2 ] + e2@eqLogLik[ e2@reporters, columns2 ]
  
  LocalLogLikelihoods( reporters = reporters
                     , actors = actors
                     , difLogLik = difLogLik
                     , eqLogLik = eqLogLik )  
})


#' Make a LocalLogLikelihoods object from a limma fit
#' 
#' @param fit - the fit object produced by limma
#' @param actors - the list of actor names
#' @param wt - string specifying wild type
#' @return a LocalLogLikelihoods object
makeLLL = function ( fit, actors, wt = "WT" ){
  
  # Extract some useful info
  theColnames = colnames(fit$coefficients)
  doubleSpecs = NULL
    for( gene1 in actors )
      for( gene2 in actors )
        if( (str = paste0(gene1,gene2)) %in% theColnames )
          doubleSpecs = append( doubleSpecs, list( list(str, c(gene1,gene2)) ) )
  
  #nActors = length( actors );
  reporters = names( fit$Amean );
  #nReporters = length( reporters );
  
  # We'll only keep reporters w/ significant changes in the WT, get them here.
  sigReporters = ( 0 != decideTests(
    eBayes( fit )[,which( theColnames == wt )] ) )
  
  # We make one contrast for each singke KO and two contrasts for each double KO.
  contrastMatrix = mkContrastMatrix( actors, doubleSpecs, wt, theColnames )
  
  # We used to... # Treat a KO effect that is opposite to the WT effect as impossible.
  # directionMismatch =
  #   sign( fit$coefficients[,wt] ) != sign( fit$coefficients[,actors] )
  
  # Get the fit model to get the log likelihoods
  fit2 = eBayes( contrasts.fit( fit, contrastMatrix ) )
  lls = getLogLikelihoods( fit2 )
  
  # return
  LocalLogLikelihoods( reporters = reporters
                     , actors = actors
                     , eqLogLik = lls$equal
                     , difLogLik = lls$different )
}

#' Function for getting log-likelihoods out of the fit
#' 
#' @param fit - limma's fit object
#' @return a list of two matrices representing the log likelihoods of equal and differential expression for each
#' reporter for each contrast in the fit.
getLogLikelihoods = function( fit ){
  
  # Compute likelihood given (t=0) from the t pdf
  lleq = dt( fit$t, fit$df.total, log = TRUE )
  
  # Compute likelihood given (t!=0) be de-posteriorifying log-odds
  lldiff = fit$lods + log(1/fit$proportion - 1) + lleq
  
  # return
  list( equal = lleq, different = lldiff )
}

#' Function for accessing a double-KO contrast in a table
#' 
#' @param table - the table to access
#' @param double1 - an actor in the double KO
#' @param double2 - another actor in the double KO
#' @param single - the singke KO/WT identifier
#' @return the column in the table matching the description, NA if not found
accessTable = function ( table, double1, double2, single ) {
  if ( ( ( refstr = paste0( double1, double2, '-', single ) ) %in% colnames( table ) ) ||
       ( ( refstr = paste0( double2, double1, '-', single ) ) %in% colnames( table ) ) )
    table[, refstr ]
  else
    NA
}