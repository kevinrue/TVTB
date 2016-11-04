
# UCSC text input is first parsed to DataFrame
# and must be checked before making a GRanges

validateDataFrameGRanges <- function(rawData, EnsDbPkg){
    validate(need(
        all(dim(rawData) >= c(1, 3)),
        "BED file must have at least 1 row and 3 columns"))

    validate(
        need(
            all(
                rawData[,1] %in%
                    seqlevels(EnsDbPkg)),
            "First column contains invalid chromosome names"),
        need(
            is.numeric(rawData[,2]),
            "Second column is not numeric"),
        need(
            is.numeric(rawData[,3]),
            "Third column is not numeric")
        )

}
