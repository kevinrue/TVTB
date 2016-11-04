
# Constructors ----

setMethod(
    f = "initialize",
    signature = c("VcfFilterList"),
    definition = function(
        .Object, exprs = list(), ...,
        active = rep(TRUE, length(c(list(...), exprs)))){

        # Fill slots with data
        .Object@filterRules <- c(list(...), exprs)
        .Object@active <- active

        validObject(.Object)

        return(.Object)
    })


setMethod(
    f = "VcfFilterList",

    definition = function(
        exprs = list(), ...,
        active = rep(TRUE, length(c(list(...), exprs)))){
        new(Class = "VcfFilterList", c(exprs, ...), active = active)
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

# Subsetting ----

setMethod(
    "[", c("VcfFilterList", "ANY", "missing", "ANY"),
    function(x, i, ...)
    {
        VcfFilterList(exprs = slot(x, "filterRules")[i])
    }
)

setMethod(
    "[[", c("VcfFilterList", "ANY", "missing"),
    function(x, i, ...)
    {
        VcfFilterList(exprs = slot(x, "filterRules")[[i]])
    }
)

# Other methods ----

### activate
# setMethod(
#     f = "activate",
#     signature = c("VcfFilterList", "numeric"),
#     definition = function(x)
#         slot(x, "active")
# )
