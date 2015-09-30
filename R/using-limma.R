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
  
  # Useful
  the.colnames = colnames(fit$coefficients);
  
  # We make one contrast for each singke KO and two contrasts for each double KO.
  contrast.matrix = mk.contrast.matrix( single.genes, double.specs, wt.str, the.colnames )
  
  log.probs$fit = eBayes( contrasts.fit(fit,contrast.matrix) );
  
  # Extract log-odds from fit.
  
  # Return.
  log.probs
}

mk.contrast.matrix = function( single.genes, double.specs, wt.str, the.colnames ){
  
  n.1le = length(single.genes);
  n.2le = length(double.specs);
  
  # Contrast generation:
  # We make one contrast for each singke KO and two contrasts for each double KO.
  contrast.matrix = matrix(data = 0, nrow = length(the.colnames), ncol = n.1le+2*n.2le );
  # Single KO contrasts:
  i = 0
  for( gene in single.genes ){
    i = i+1
    contrast.matrix[which(the.colnames == gene),i] = 1;
    contrast.matrix[which(the.colnames == wt.str),i] = -1;
  }
  # Double KO contrasts:
  for( element in double.specs )
    for( single in element[[2]] ){
      i = i+1
      contrast.matrix[which(the.colnames == element[[1]]),i] = 1;
      contrast.matrix[which(the.colnames == single),i] = -1;
    }
  
  contrast.matrix
}