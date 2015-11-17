# Class definition
LimmaLogProbs = setClass("LimmaLogProbs",
                         slots = c(
                           actors = "character",
                           reporters = "character",
                           nActors = "numeric",
                           nReporters = "numeric",
                           prior = "numeric",
                           singleGtWT = "matrix",
                           doubleVsingle = "list"
                         ))

# Here we associate the generics from laps-interface with the class.
setMethod(f = "scoreIndependentPathways",
          signature(theObject = "LimmaLogProbs", actor1 = "character", actor2 = "character"),
          definition = function(theObject, actor1, actor2){
            (
              if( actor1 %in% names(theObject@doubleVsingle) &&
                  actor2 %in% names( theObject@doubleVsingle[[actor1]] ) )
                theObject@doubleVsingle[[actor1]][[actor2]][["gt"]]
              else -log(2)
            ) + (
              if( actor2 %in% names(theObject@doubleVsingle) &&
                  actor1 %in% names( theObject@doubleVsingle[[actor2]] ) )
                theObject@doubleVsingle[[actor2]][[actor1]][["gt"]]
              else -log(2)
            )
          })

setMethod(f = "scoreIndependentPathways",
          signature(theObject = "LimmaLogProbs", actor1 = "numeric", actor2 = "numeric"),
          definition = function(theObject,actor1,actor2){
            scoreIndependentPathways(theObject, theObject@actors[actor1],theObject@actors[actor2])
          })

setMethod(f = "scoreSharedPathways",
          signature(theObject = "LimmaLogProbs", actor1 = "character", actor2 = "character"),
          definition = function(theObject, actor1, actor2){
            (
              if( actor1 %in% names(theObject@doubleVsingle) &&
                  actor2 %in% names( theObject@doubleVsingle[[actor1]] ) )
                theObject@doubleVsingle[[actor1]][[actor2]][["eq"]]
              else -log(2)
            ) + (
              if( actor2 %in% names(theObject@doubleVsingle) &&
                  actor1 %in% names( theObject@doubleVsingle[[actor2]] ) )
                theObject@doubleVsingle[[actor2]][[actor1]][["eq"]]
              else -log(2)
            )
          })

setMethod(f = "scoreSharedPathways",
          signature(theObject = "LimmaLogProbs", actor1 = "numeric", actor2 = "numeric"),
          definition = function(theObject,actor1,actor2){
            scoreSharedPathways(theObject, theObject@actors[actor1],theObject@actors[actor2])
          })

setMethod(f = "ancestryScoreMatrix",
          signature = "LimmaLogProbs",
          definition = function(theObject) t(theObject@singleGtWT) )

setMethod(f = "getActors",
          signature = "LimmaLogProbs",
          definition = function(theObject) theObject@actors )

setMethod(f = "getReporters",
          signature = "LimmaLogProbs",
          definition = function(theObject) theObject@reporters )

setMethod(f= "howManyActors",
          signature = "LimmaLogProbs",
          definition = function(theObject) theObject@nActors )

setMethod(f= "howManyReporters",
          signature = "LimmaLogProbs",
          definition = function(theObject) theObject@nReporters )

# Use limma to get lprobs.
# Starting point is having a 'fit' (before contrasts).
# Automatically build contrasts and run the right models.
makeLimmaLogProbs = function(
  fit,
  actors,
  theColnames = colnames(fit$coefficients),
  doubleSpecs = {
    ds = NULL
    for( gene1 in actors )
      for( gene2 in actors )
        if( (str = paste0(gene1,gene2)) %in% theColnames )
          ds = append( ds, list( list(str, c(gene1,gene2)) ) )
    ds
  },
  wt = "WT", prior = 0.01 ){
  # fit is the model fit object from limma obtained with lmFit
  # single.genes is a list of strings identifying the single genes
  # double.specs is a list where each element has a list with a first element identifying the name of a double KO,
  #             and the second being a pair of the corresponding single KOs. E.g. "hog1msn2", ["hog1", "msn2"]
  
  nActors = length( actors );
  reporters = names(fit$Amean);
  nReporters = length( reporters );
  
  # We'll only keep reporters w/ significant changes in the WT, get them here.
  sigReporters = ( 0 != decideTests(
    eBayes( fit, proportion = prior )[,which( theColnames == wt )] ) )

  # We make one contrast for each singke KO and two contrasts for each double KO.
  contrastMatrix = mkContrastMatrix( actors, doubleSpecs, wt, theColnames )
  
  # Treat a KO effect that is opposite to the WT effect as impossible.
  directionMismatch =
    sign( fit$coefficients[,wt] ) != sign( fit$coefficients[,actors] )
  
  # Get the fit model to get the lods matrix
  lodsMatrix = eBayes( contrasts.fit(fit,contrastMatrix), proportion = prior )$lods
  
  # Extract log-probs from fit.
  singleGtWT = lprob.from.lods( lodsMatrix[,1:nActors] )
  singleGtWT[directionMismatch] = -Inf
  
  # We don't want to deal with NA "single>WT" probabilities, so those reporters that
  # yield those will be thrown out too.
  sigReporters = sigReporters & !apply(is.na(singleGtWT),1,any)
  
  # Now we trim
  singleGtWT = singleGtWT[sigReporters,]
  
  doubleVsingle = NULL

  i = nActors
  for( element in doubleSpecs ){
    gene = element[[2]]
    for( j in 1:2 ){
      i = i+1
      
      lprobs = lprob.from.lods( lodsMatrix[sigReporters,i] )
      
      directionMismatch =
        sign( fit$coefficients[sigReporters,gene[j]] ) !=
        sign( lodsMatrix[sigReporters,i] )
      
      item=NULL
      item$eq = log1mexp( lprobs )
      item$gt = lprobs
      item$gt[directionMismatch] = -Inf
      doubleVsingle[[gene[j]]][[gene[3-j]]] = item
    }
  }
  
  # Update reporter list
  reporters = reporters[sigReporters]
  nReporters = length(reporters)
  
  # Return.
  LimmaLogProbs(
    actors = actors,
    nActors = nActors,
    reporters = reporters,
    nReporters = nReporters,
    prior = prior,
    singleGtWT = singleGtWT,
    doubleVsingle = doubleVsingle
  )
}

mkContrastMatrix = function( single.genes, double.specs, wt.str, the.colnames ){
  
  n.1le = length(single.genes)
  n.2le = length(double.specs)
  
  # Contrast generation:
  # We make one contrast for each singke KO and two contrasts for each double KO.
  contrast.matrix = matrix(data = 0, nrow = length(the.colnames), ncol = n.1le+2*n.2le )
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
    }
  
  contrast.matrix
}

setMethod("+", signature(e1 = "LimmaLogProbs", e2 = "LimmaLogProbs"), function (e1,e2){
  
  if( e1@prior != e2@prior ) stop("Mismatched priors when merging limma log-prob objects")
  
  reporters = intersect( e1@reporters, e2@reporters )
  r1s = which( e1@reporters %in% reporters )
  r2s = which( e2@reporters %in% reporters )
  nReporters = length(reporters)
  prior = e1@prior
  
  actors.intersection = intersect( e1@actors, e2@actors )
  actors.only.1 = setdiff( e1@actors, actors.intersection )
  actors.only.2 = setdiff( e2@actors, actors.intersection )

  actors = c( actors.only.1, actors.only.2, actors.intersection )

  nActors = length( actors )
  
  # Merge the single>WT matrix
  singleGtWT = cbind(
    e1@singleGtWT[r1s,which(e1@actors %in% actors.only.1)],
    e2@singleGtWT[r2s,which(e2@actors %in% actors.only.2)],
    combinePosteriors(
      e1@singleGtWT[r1s,which(e1@actors %in% actors.intersection)],
      e2@singleGtWT[r2s,which(e2@actors %in% actors.intersection)],
      prior ) )

  # Merge the double vs single matrix
  doubleVsingle = e1@doubleVsingle
  
  # Clear reporters
  for ( l1 in e1@doubleVsingle )
    for ( l2 in l1 ){
      l2$gt = l2$gt[r1s]
      l2$eq = l2$eq[r1s]
    }
  
  for ( gene1 in names(e2@doubleVsingle) )
    for ( gene2 in names(e2@doubleVsingle[[gene1]]) ){
      if ( is.null( doubleVsingle[[gene1]][[gene2]] ) ){
        doubleVsingle[[gene1]][[gene2]]$gt = e2@doubleVsingle[[gene1]][[gene2]]$gt[r2s]
        doubleVsingle[[gene1]][[gene2]]$eq = e2@doubleVsingle[[gene1]][[gene2]]$eq[r2s]
      } else {
        doubleVsingle[[gene1]][[gene2]]$gt = combinePosteriors(
          e2@doubleVsingle[[gene1]][[gene2]]$gt[r2s],
          doubleVsingle[[gene1]][[gene2]]$gt,
          prior )
        doubleVsingle[[gene1]][[gene2]]$eq = combinePosteriors(
          e2@doubleVsingle[[gene1]][[gene2]]$eq[r2s],
          doubleVsingle[[gene1]][[gene2]]$eq,
          1-prior )
      }
    }
  
  LimmaLogProbs(
    actors = actors,
    reporters = reporters,
    prior = prior,
    nActors = nActors,
    nReporters = nReporters,
    singleGtWT = singleGtWT,
    doubleVsingle = doubleVsingle
  )
})

combinePosteriors = function( p1, p2, prior ){
  lprob.from.lods( lprob2lods( p1 ) + lprob2lods( p2 ) - log(prior/(1-prior)) )
}