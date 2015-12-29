# laps-interface.R
# Defines the R4 "interface" for the local ancestry probability scores
# I'm guessing in R the interface is just generics.
#' @export scoreIndependentPathways scoreSharedPathways ancestryScoreMatrix getActors getReporters howManyActors howManyReporters nonAncestryScoreMatrix scoreSingleAncestor scoreNeitherAncestor

setGeneric(name = "scoreIndependentPathways",
           def = function(theObject, actor1, actor2)
           { standardGeneric("scoreIndependentPathways") } )
setGeneric(name = "scoreSharedPathways",
           def = function(theObject, actor1, actor2)
           { standardGeneric("scoreSharedPathways") } )
setGeneric(name = "ancestryScoreMatrix",
           def = function(theObject)
           { standardGeneric("ancestryScoreMatrix") } )
setGeneric(name = "getActors",
           def = function(theObject)
           { standardGeneric("getActors") } )
setGeneric(name = "getReporters",
           def = function(theObject)
           { standardGeneric("getReporters") } )
setGeneric(name = "howManyActors",
           def = function(theObject)
           { length( getActors( theObject ) ) } )
setGeneric(name = "howManyReporters",
           def = function(theObject)
           { length( getReporters( theObject ) ) } )

# These were added with the new likelihood-based object
setGeneric( name = "nonAncestryScoreMatrix"
          , def = function ( theObject ) {
            standardGeneric( "nonAncestryScoreMatrix")
          } )
setGeneric( name = "scoreSingleAncestor"
          , def = function ( theObject, ancestor, nonancestor ) {
            standardGeneric( "scoreSingleAncestor" )
          } )
setGeneric( name = "scoreNeitherAncestor"
            , def = function ( theObject, actor1, actor2 ) {
              standardGeneric( "scoreNeitherAncestor" )
            } )