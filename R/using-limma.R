# Use limma to get lprobs.
# Starting point is having a 'fit' (before contrasts).
# Automatically build contrasts and run the right models.
log.probs = function( fit, single.genes, double.specs, wt.str = "WT" ){
  # fit is the model fit object from limma obtained with lmFit
  # single.genes is a list of strings identifying the single genes
  # double.specs is a list where each element has a list with a first element identifying the name of a double KO,
  #             and the second being a pair of the corresponding single KOs. E.g. "hog1msn2", ["hog1", "msn2"]
  
  log.probs = NULL;
  log.probs$actors = single.genes;
  log.probs$n.actors = length( single.genes );
  log.probs$reporters = names(fit$Amean);
  log.probs$n.reporters = length( log.probs$reporters );
  
  # We make one contrast for each singke KO and two contrasts for each double KO.
  contrast.matrix = mk.contrast.matrix( single.genes, double.specs, wt.str, the.colnames )
  
  direction.mismatch =
    sign( fit$coefficients[,wt.str] ) != sign( fit$coefficients[,single.genes] )
    
  fit.for.contrasts = eBayes( contrasts.fit(fit,contrast.matrix) );
  
  # Extract log-probs from fit.
  log.probs$single.gt.wt = lprob.from.lods( fit.for.contrasts$lods[,1:log.probs$n.actors] )
  log.probs$single.gt.wt[direction.mismatch] = -Inf
  
  log.probs$double.vs.single = NULL

  i = log.probs$n.actors
  for( element in double.specs ){
    gene = element[[2]]
    for( j in 1:2 ){
      i = i+1
      
      lprobs = lprob.from.lods( fit.for.contrasts$lods[,i] )
      
      direction.mismatch =
        sign( fit$coefficients[,gene[j]] ) !=
        sign( fit.for.contrasts$coefficients[,i] )
      
      item=NULL
      item$eq = log1mexp( lprobs )
      item$gt = lprobs
      item$gt[direction.mismatch] = -Inf
      log.probs$double.vs.single[[gene[j]]][[gene[3-j]]] = item
    }
  }
  
  # Return.
  log.probs
}

mk.contrast.matrix = function( single.genes, double.specs, wt.str, the.colnames ){
  
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