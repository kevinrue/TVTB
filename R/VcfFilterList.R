
# Constructors ----

setMethod(
    f = "initialize",
    signature = c("VcfFilterList"),
    definition = function(
        .Object, ...){

        # Fill slots with data
        .Object@filterRules <- list(...)

        validObject(.Object)

        return(.Object)
    })


setMethod(
    f = "VcfFilterList",

    definition = function(...){

        new(Class = "VcfFilterList", ...)
    }
)

# Getters and Setters ----

### filterRules
setMethod(
    f = "filterRules",
    signature = c("VcfFilterList"),
    definition = function(x)
        slot(x, "filterRules")
)

setReplaceMethod(
    f = "filterRules", c("VcfFilterList", "VcfFilterList"),
    function(x, value){
        slot(x, "filterRules") <- value
        validObject(x)
        return(x)
    }
)

setReplaceMethod(
    f = "filterRules", c("VcfFilterList", "VcfFixedFilter"),
    function(x, value){
        slot(x, "filterRules") <- list(value)
        validObject(x)
        return(x)
    }
)

setReplaceMethod(
    f = "filterRules", c("VcfFilterList", "VcfInfoFilter"),
    function(x, value){
        slot(x, "filterRules") <- list(value)
        validObject(x)
        return(x)
    }
)

setReplaceMethod(
    f = "filterRules", c("VcfFilterList", "VcfInfoFilter"),
    function(x, value){
        slot(x, "filterRules") <- list(value)
        validObject(x)
        return(x)
    }
)
