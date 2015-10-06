# Class definition
LimmaLogProbs = setClass("LimmaLogProbs",
                         slots = c(
                           actors = "character",
                           reporters = "character",
                           nActors = "numeric",
                           nReporters = "numeric",
                           singleGtWT = "numeric",
                           doubleVsingle = "list"
                         ))

# Use limma to get lprobs.
# Starting point is having a 'fit' (before contrasts).
# Automatically build contrasts and run the right models.
llp = function( fit, actors, doubleSpecs, wt = "WT" ){
  # fit is the model fit object from limma obtained with lmFit
  # single.genes is a list of strings identifying the single genes
  # double.specs is a list where each element has a list with a first element identifying the name of a double KO,
  #             and the second being a pair of the corresponding single KOs. E.g. "hog1msn2", ["hog1", "msn2"]
  
  nActors = length( actors );
  reporters = names(fit$Amean);
  nReporters = length( reporters );
  theColnames = colnames(fit$coefficients);
  
  # We make one contrast for each singke KO and two contrasts for each double KO.
  contrastMatrix = mkContrastMatrix( actors, doubleSpecs, wt, theColnames )
  
  directionMismatch =
    sign( fit$coefficients[,wt] ) != sign( fit$coefficients[,actors] )
    
  fit4contrasts = eBayes( contrasts.fit(fit,contrastMatrix) );
  
  # Extract log-probs from fit.
  singleGtWT = lprob.from.lods( fit4contrasts$lods[,1:nActors] )
  singleGtWT[directionMismatch] = -Inf
  
  doubleVsingle = NULL

  i = nActors
  for( element in doubleSpecs ){
    gene = element[[2]]
    for( j in 1:2 ){
      i = i+1
      
      lprobs = lprob.from.lods( fit4contrasts$lods[,i] )
      
      directionMismatch =
        sign( fit$coefficients[,gene[j]] ) !=
        sign( fit4contrasts$coefficients[,i] )
      
      item=NULL
      item$eq = log1mexp( lprobs )
      item$gt = lprobs
      item$gt[directionMismatch] = -Inf
      log.probs$doubleVsingle[[gene[j]]][[gene[3-j]]] = item
    }
  }
  
  # Return.
  LimmaLogProbs(
    actors = actors,
    nActors = nActors,
    reporters = reporters,
    nReporters = nReporters,
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

merge.log.probs = function( lp1, lp2, prior = 0.01 ){
  result = NULL
  
  if( any( lp1$reporters != lp2$reporters) ) stop("Mismatching reporter lists not yet supported")
  
  result$reporters = lp1$reporters
  result$n.reporters = lp1$n.reporters
  
  actors.intersection = intersect( lp1$actors, lp2$actors )
  actors.only.1 = setdiff( lp1$actors, actors.intersection )
  actors.only.2 = setdiff( lp2$actors, actors.intersection )

  result$actors = c( actors.only.1, actors.only.2, actors.intersection )

  result$n.actors = length( result$actors )
  
  result$single.gt.wt = cbind( lp1$single.gt.wt[,which(lp1$actors %in% actors.only.1)],
                               lp2$single.gt.wt[,which(lp2$actors %in% actors.only.2)],
                               lp1$single.gt.wt[,which(lp1$actors %in% actors.intersection)] +
                                 lp2$single.gt.wt[,which(lp2$actors %in% actors.intersection)] -
                                 log(prior/(1-prior)))

  result$double.vs.single = lp1$double.vs.single
  
  for( gene1 in names(lp2$double.vs.single) )
    for( gene2 in names(lp2$double.vs.single[[gene1]])){
      if( is.null( result$double.vs.single[[gene1]][[gene2]] ))
        result$double.vs.single[[gene1]][[gene2]] = lp2$double.vs.single[[gene1]][[gene2]]
      else
        result$double.vs.single[[gene1]][[gene2]] = lp2$double.vs.single[[gene1]][[gene2]] + result$double.vs.single[[gene1]][[gene2]] - log(prior/(1-prior))
    }
  
  result
}