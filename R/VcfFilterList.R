
# Constructors ----

setMethod(
    f = "initialize",
    signature = c("VcfFilterList"),
    definition = function(
        .Object, ..., active = rep(TRUE, length(list(...)))){

        # Fill slots with data
        .Object@filterRules <- list(...)
        .Object@active <- active

        validObject(.Object)

        return(.Object)
    })


setMethod(
    f = "VcfFilterList",

    definition = function(..., active = rep(TRUE, length(list(...)))){
        new(Class = "VcfFilterList", ..., active = active)
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

### active
setMethod(
    f = "active",
    signature = c("VcfFilterList"),
    definition = function(x)
        slot(x, "active")
)

setReplaceMethod(
    f = "active", c("VcfFilterList"),
    function(x, value){

        if (is.logical(value)){
            slot(x, "active") <- value
        } else {
            stop("Invalid value class: ", class(value))
        }

        validObject(x)
        return(x)
    }
)
