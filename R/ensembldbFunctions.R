
## Based on the given EnsDb package name it loads the library and returns
## the object.
getEdb <- function(x){
    stopifnot(requireNamespace(x))
    return(get(x = x, envir = asNamespace(x)))
}

# If "," present, taken as the value separator and " " are trimmed
# Otherwise, " " taken as the value separator
# Returns an ensembldb::SeqendFilter of the appropriate type
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
