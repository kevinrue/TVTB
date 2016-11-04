
tryParseCsq <- function(vcf, csqField){

    message("Parsing ", csqField," to GRanges ...")

    rawCsq <- tryCatch(
        parseCSQToGRanges(
            x = vcf,
            VCFRowID = rownames(vcf),
            info.key = csqField),
        warning = function(warn){
            warning(warn)
            return(NULL)},
        error = function(err){
            warning(geterrmessage())
            return(NULL)
        }
    )
}

tryParseBed <- function(bed){
    message("Parsing ", bed," as GRanges ...")

    rawCsq <- tryCatch(
        import.bed(con = bed),
        warning = function(warn){
            warning(warn)
            return(NULL)},
        error = function(err){
            warning(geterrmessage())
            return(NULL)
        }
    )
}
