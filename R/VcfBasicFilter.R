
# Constructors ----

setMethod(
    f = "initialize",
    signature = c("VcfFixedFilter"),
    definition = function(
        .Object, name, condition, value){

        # Fill slots with data
        .Object@name <- name
        .Object@condition <- condition
        .Object@value <- value

        validObject(.Object)

        return(.Object)
    })


setMethod(
    f = "VcfFixedFilter",
    signature = c(name="character", condition="character", value="ANY"),
    definition = function(
        name, condition, value){

        new(
            Class = "VcfFixedFilter",
            name = name, condition = condition, value = value)
    }
)

setMethod(
    f = "initialize",
    signature = c("VcfInfoFilter"),
    definition = function(
        .Object, name, condition, value){

        # Fill slots with data
        .Object@name <- name
        .Object@condition <- condition
        .Object@value <- value

        validObject(.Object)

        return(.Object)
    })


setMethod(
    f = "VcfInfoFilter",
    signature = c(name="character", condition="character", value="ANY"),
    definition = function(
        name, condition, value){

        new(
            Class = "VcfInfoFilter",
            name = name, condition = condition, value = value)
    }
)

setMethod(
    f = "initialize",
    signature = c("VcfVepFilter"),
    definition = function(
        .Object, name, condition, value){

        # Fill slots with data
        .Object@name <- name
        .Object@condition <- condition
        .Object@value <- value

        validObject(.Object)

        return(.Object)
    })


setMethod(
    f = "VcfVepFilter",
    signature = c(name="character", condition="character", value="ANY"),
    definition = function(
        name, condition, value){

        new(
            Class = "VcfVepFilter",
            name = name, condition = condition, value = value)
    }
)

# Getters and Setters ----

### name
setMethod(
    f = "name",
    signature = c("VcfFixedFilter"),
    definition = function(x)
        slot(x, "name")
)

setMethod(
    f = "name",
    signature = c("VcfInfoFilter"),
    definition = function(x)
        slot(x, "name")
)

setMethod(
    f = "name",
    signature = c("VcfVepFilter"),
    definition = function(x)
        slot(x, "name")
)

setReplaceMethod(
    f = "name", c("VcfFixedFilter", "character"),
    function(x, value){
        slot(x, "name") <- value
        validObject(x)
        return(x)
    }
)

setReplaceMethod(
    f = "name", c("VcfInfoFilter", "character"),
    function(x, value){
        slot(x, "name") <- value
        validObject(x)
        return(x)
    }
)

setReplaceMethod(
    f = "name", c("VcfVepFilter", "character"),
    function(x, value){
        slot(x, "name") <- value
        validObject(x)
        return(x)
    }
)

### condition
setMethod(
    f = "condition",
    signature = c("VcfFixedFilter"),
    definition = function(x)
        slot(x, "condition")
)

setMethod(
    f = "condition",
    signature = c("VcfInfoFilter"),
    definition = function(x)
        slot(x, "condition")
)

setMethod(
    f = "condition",
    signature = c("VcfVepFilter"),
    definition = function(x)
        slot(x, "condition")
)

setReplaceMethod(
    f = "condition", c("VcfFixedFilter", "character"),
    function(x, value){
        slot(x, "condition") <- value
        validObject(x)
        return(x)
    }
)

setReplaceMethod(
    f = "condition", c("VcfInfoFilter", "character"),
    function(x, value){
        slot(x, "condition") <- value
        validObject(x)
        return(x)
    }
)

setReplaceMethod(
    f = "condition", c("VcfVepFilter", "character"),
    function(x, value){
        slot(x, "condition") <- value
        validObject(x)
        return(x)
    }
)

### value
setMethod(
    f = "value",
    signature = c("VcfFixedFilter"),
    definition = function(x)
        slot(x, "value")
)

setMethod(
    f = "value",
    signature = c("VcfInfoFilter"),
    definition = function(x)
        slot(x, "value")
)

setMethod(
    f = "value",
    signature = c("VcfVepFilter"),
    definition = function(x)
        slot(x, "value")
)

setReplaceMethod(
    f = "value", c("VcfFixedFilter", "ANY"),
    function(x, value){
        slot(x, "value") <- value
        validObject(x)
        return(x)
    }
)

setReplaceMethod(
    f = "value", c("VcfInfoFilter", "ANY"),
    function(x, value){
        slot(x, "value") <- value
        validObject(x)
        return(x)
    }
)

setReplaceMethod(
    f = "value", c("VcfVepFilter", "ANY"),
    function(x, value){
        slot(x, "value") <- value
        validObject(x)
        return(x)
    }
)
