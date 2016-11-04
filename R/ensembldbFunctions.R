
## Based on the given EnsDb package name it loads the library and returns
## the object.
getEdb <- function(x){
    require(x, character.only = TRUE)
    return(get(x))
}

EnsDbFilter <- function(type, condition, value){
    ## split by "," if present
    if(length(grep(x = value, pattern = ",")) > 0){
        ##  remove whitespaces if "," are present
        value <- gsub(value, pattern = " ", replacement = "", fixed = TRUE)
        value <- unlist(strsplit(value, split = ","))
    }
    if(length(grep(x = value, pattern = " ", fixed = TRUE)) > 0){
        value <- unlist(strsplit(x = value, split = " ", fixed = TRUE))
    }

    filter <- switch(
        type,
        Genename = GenenameFilter(value = value, condition = condition)
    )

    if (is.null(filter)){
        warning("Invalid filter type")
    }

    return(filter)
}
