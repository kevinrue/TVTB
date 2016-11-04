
.infoMinFilter <- function(vcf, key, value){
    info(vcf)[,key] >= value
}

.infoMaxFilter <- function(vcf, key, value){
    info(vcf)[,key] <= value
}

.applyVcfFilters <- function(
    vcf,
    MAF = NULL){

    df <- DataFrame(row.names = rownames(vcf))

    # That could be parallelised
    if (!is.null(MAF)){
        stopifnot(is.numeric(MAF))
        stopifnot(length(MAF) == 2)

        stopifnot(MAF[1] >= 0)
        stopifnot(MAF[2] <= 0.5)
        stopifnot(MAF[2] >= MAF[1])

        df[,"mafMin"] <- .infoMinFilter(vcf = vcf, key = "MAF", value = MAF[1])
        df[,"mafMax"] <- .infoMaxFilter(vcf = vcf, key = "MAF", value = MAF[2])
    }

    return(df)
}

.extractFilterRules <- function(filterDataFrame){
    FilterRules(
        sapply(
            X = colnames(filterDataFrame),
            FUN = function(x){parse(text = x)}, simplify = FALSE
        )
    )
}

filterExpandedVcf <- function(
    vcf,
    MAF = NULL){
    message("filterExpandedVcf...")

    # Boolean DataFrame of boolean values answering tests
    filterDataFrame <- .applyVcfFilters(
        vcf = vcf,
        MAF = MAF)

    # One column for each test
    if (ncol(filterDataFrame) > 0){
        # Make a FilterRules object from column names
        filterRules <- .extractFilterRules(filterDataFrame = filterDataFrame)
        # Evaluate filters (equivalent to rowSums = ncol)
        # TODO: return DataFrame and let Shiny eval the filters with
        # checkboxes turning filters on/off.
        filterResults <- S4Vectors::eval(
            expr = filterRules,
            envir = filterDataFrame)
        # Return VCF records passing filters
        return(vcf[which(filterResults)])
    } else {
        return(vcf)
    }

}
