# Instantiate a BasicFilter of the appropriate class
EnsDbFilter <- function(type, condition, value){

    # Use "," as separator if present
    if(length(grep(",", value)) > 0){
        # Remove whitespaces if "," are used as separator
        value <- gsub(" ", "", value, fixed = TRUE)
        value <- unlist(strsplit(value, split = ","))
    }

    # Use " " as separator if they were not trimmed above (no "," present)
    if(length(grep(x = value, pattern = " ", fixed = TRUE)) > 0){
        value <- unlist(strsplit(value, " ", fixed = TRUE))
    }

    # Instantiate a filter of the appropriate class
    # TODO: support more classes of filters
    filter <- switch(
        type,
        Genename = ensembldb::GenenameFilter(value, condition),
        stop("Invalid filter type")
    )

    return(filter)
}
