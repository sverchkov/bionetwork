# This is the new likelihood-based scoring class and its functions
library( limma )

# Class definition
LocalLogLikelihoods = setClass( "LocalLogLikelihoods",
                                slots = c( actors = "character"
                                         , reporters = "character"
                                         , eqLogLik = "matrix" # Log-likelihoods given equal expression change
                                         , difLogLik = "matrix" # Log-likelihoods given different expression change
                                        #, ampLogLik = "matrix" # Log-likelihoods given amplified expression change
                                ) )

# Simple property functions
setMethod(f = "getActors",
          signature = "LocalLogLikelihoods",
          definition = function(theObject) theObject@actors )

setMethod(f = "getReporters",
          signature = "LocalLogLikelihoods",
          definition = function ( theObject ) theObject@reporters )

# Functions to get log-likelihoods under equal expression
setMethod(f = "getLLGivenEqual",
          signature( theObject = "LocalLogLikelihoods", actor = "character" ),
          definition = function ( theObject, actor ) theObject@eqLogLik[, actor ] )

setMethod(f = "getLLGivenDifferent",
          signature( theObject = "LocalLogLikelihoods", actor = "character" ),
          definition = function ( theObject, actor ) theObject@difLogLik[, actor ] )

setMethod(f = "getLLGivenEqual2",
          signature( theObject = "LocalLogLikelihoods", double1 = "character", double2 = "character", single = "character" ),
          definition = function ( theObject, double1, double2, single ){
            if ( !( ( str = paste0( double1, double2, '-', single ) ) %in% colnames( eqLogLik ) ) &&
                 !( ( str = paste0( double2, double1, '-', single ) ) %in% colnames( eqLogLik ) ) )
              0
            else
              theObject@eqLogLik[, str ]
          } )

setMethod(f = "getLLGivenDifferent2",
          signature( theObject = "LocalLogLikelihoods", double1 = "character", double2 = "character", single = "character" ),
          definition = function ( theObject, double1, double2, single ){
            if ( !( ( str = paste0( double1, double2, '-', single ) ) %in% colnames( eqLogLik ) ) &&
                 !( ( str = paste0( double2, double1, '-', single ) ) %in% colnames( eqLogLik ) ) )
              0
            else
              theObject@difLogLik[, str ]
          } )

# We might write in numeric equivalents later. For now the char enforcement should do.

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
    eBayes( fit, proportion = prior )[,which( theColnames == wt )] ) )
  
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