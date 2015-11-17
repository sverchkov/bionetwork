# laps-interface.R
# Defines the R4 "interface" for the local ancestry probability scores
# I'm guessing in R the interface is just generics.

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
           { standardGeneric("howManyActors") } )
setGeneric(name = "howManyReporters",
           def = function(theObject)
           { standardGeneric("howManyReporters") } )
