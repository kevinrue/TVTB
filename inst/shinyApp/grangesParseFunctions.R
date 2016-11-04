
# UCSC text input is first parsed to DataFrame
# and must be checked before making a GRanges

validateDataFrameGRanges <- function(rawData, EnsDbPkg){

  if (is.null(EnsDbPkg)){
    stop("An annotation package must be supplied to validate seqnames.")
  }

  # base::dim called, as rawData is data.frame
  if (!all(dim(rawData) >= c(1, 3))){
    stop("BED file must have at least 1 row and 3 columns")
  }

  if (!all(rawData[,1] %in% ensembldb::seqlevels(EnsDbPkg))){
    stop("First column contains invalid chromosome names")
  }

  if (!is.numeric(rawData[,2])){
    stop("Second column is not numeric")
  }

  if (!is.numeric(rawData[,3])){
    stop("Third column is not numeric")
  }

  if (!all(rawData[,3] >= rawData[,2])){
    stop("All ranges must have positive length")
  }

  return(GenomicRanges::GRanges(
    seqnames = rawData[,1],
    ranges = IRanges::IRanges(
      start = rawData[,2],
      end = rawData[,3]
    )
  ))
}
