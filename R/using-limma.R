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

  # We make one contrast for each singke KO and two contrasts for each double KO.
  contrastMatrix = mkContrastMatrix( actors, doubleSpecs, wt, theColnames )
  
  directionMismatch =
    sign( fit$coefficients[,wt] ) != sign( fit$coefficients[,actors] )
    
  fit4contrasts = eBayes( contrasts.fit(fit,contrastMatrix), proportion = prior );
  
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
      doubleVsingle[[gene[j]]][[gene[3-j]]] = item
    }
  }
  
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
  
  if( any( e1@reporters != e2@reporters) ) stop("Mismatching reporter lists not yet supported")
  if( e1@prior != e2@prior ) stop("Mismatched priors when merging limma log-prob objects")
  
  reporters = e1@reporters
  nReporters = e1@nReporters
  prior = e1@prior
  
  actors.intersection = intersect( e1@actors, e2@actors )
  actors.only.1 = setdiff( e1@actors, actors.intersection )
  actors.only.2 = setdiff( e2@actors, actors.intersection )

  actors = c( actors.only.1, actors.only.2, actors.intersection )

  nActors = length( actors )
  
  singleGtWT = cbind( e1@singleGtWT[,which(e1@actors %in% actors.only.1)],
                               e2@singleGtWT[,which(e2@actors %in% actors.only.2)],
                               e1@singleGtWT[,which(e1@actors %in% actors.intersection)] +
                                 e2@singleGtWT[,which(e2@actors %in% actors.intersection)] -
                                 log(prior/(1-prior)))

  doubleVsingle = e1@doubleVsingle
  
  for( gene1 in names(e2@doubleVsingle) )
    for( gene2 in names(e2@doubleVsingle[[gene1]])){
      if( is.null( doubleVsingle[[gene1]][[gene2]] ))
        doubleVsingle[[gene1]][[gene2]] = e2@doubleVsingle[[gene1]][[gene2]]
      else
        doubleVsingle[[gene1]][[gene2]] = e2@doubleVsingle[[gene1]][[gene2]] + doubleVsingle[[gene1]][[gene2]] - log(prior/(1-prior))
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