
# Constructors ----

setMethod(
    f = "initialize",
    signature = c("NamedFilter"),
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
    f = "NamedFilter",
    signature = c(name="character", condition="character", value="ANY"),
    definition = function(
        name, condition, value){

        new(
            Class = "NamedFilter",
            name = name, condition = condition, value = value)
    }
)

# Getters and Setters ----

### name
setMethod(
    f = "name",
    signature = c("NamedFilter"),
    definition = function(x)
        slot(x, "name")
)

setReplaceMethod(
    f = "name", c("NamedFilter", "character"),
    function(x, value){
        slot(x, "name") <- value
        validObject(x)
        return(x)
    }
)

### condition
setMethod(
    f = "condition",
    signature = c("NamedFilter"),
    definition = function(x)
        slot(x, "condition")
)

setReplaceMethod(
    f = "condition", c("NamedFilter", "character"),
    function(x, value){
        slot(x, "condition") <- value
        validObject(x)
        return(x)
    }
)

### value
setMethod(
    f = "value",
    signature = c("NamedFilter"),
    definition = function(x)
        slot(x, "value")
)

setReplaceMethod(
    f = "value", c("NamedFilter", "ANY"),
    function(x, value){
        slot(x, "value") <- value
        validObject(x)
        return(x)
    }
)
