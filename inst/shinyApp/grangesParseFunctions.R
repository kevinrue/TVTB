
# UCSC text input is first parsed to DataFrame
# and must be checked before making a GRanges

validateDataFrameGRanges <- function(rawData, EnsDbPkg){

    # base::dim called, as rawData is data.frame
    if (!all(dim(rawData) >= c(1, 3))){
        warning("BED file must have at least 1 row and 3 columns")
        return(NULL)
    }

    if (!all(rawData[,1] %in% ensembldb::seqlevels(EnsDbPkg))){
        warning("First column contains invalid chromosome names")
        return(NULL)
    }

    if (!is.numeric(rawData[,2])){
        warning("Second column is not numeric")
        return(NULL)
    }

    if (!is.numeric(rawData[,3])){
        warning("Third column is not numeric")
        return(NULL)
    }

    if (!all(rawData[,3] >= rawData[,2])){
        warning("All ranges must have positive length")
        return(NULL)
    }

    return(GenomicRanges::GRanges(
        seqnames = rawData[,1],
        ranges = IRanges::IRanges(
            start = rawData[,2],
            end = rawData[,3]
        )
    ))
}
