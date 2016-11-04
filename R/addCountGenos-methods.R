
# samples as default (all): probably the most common usage
setMethod(
    f = "addCountGenos",
    signature = c(vcf="ExpandedVCF"),
    definition = function(
        vcf, genos, key, description, samples = 1:ncol(vcf), force = FALSE){
        
        .addCountGenos(
            vcf, genos, key, description, samples = samples,
            force = force)
    
    }
)

.addCountGenos <- function(
    vcf, genos, key, description, samples = 1:ncol(vcf), force = FALSE){
    
    # Validate inputs
    stopifnot(class(samples) %in% c("integer", "numeric", "character"))
    stopifnot(is.logical(force))
    stopifnot(is.character(key))
    stopifnot(identical(x = length(key), y = as.integer(1)))
    stopifnot(is.character(description))
    stopifnot(identical(x = length(description), y = as.integer(1)))

    # Prepare the new INFO header
    newInfoHeader <- DataFrame(
        Number = 1,
        Type = "Integer",
        Description = description,
        row.names = key)

    # Identify the position of the key in the header/data if present
    keyDataIndex <- which(colnames(info(vcf)) == key)
    keyHeaderIndex <- which(rownames(info(header(vcf))) == key)

    # Check: warn/stop if key present in data
    if (length(keyDataIndex) > 0){
        if (force){
            message(key, " overwritten.")
            vcf <- .dropInfoData(vcf = vcf, key = key)
        } else {
            stop(key, " key already present in INFO.")
        }
        newData <- .countGenos(
            x = geno(vcf)[["GT"]][,samples], genos = genos)
    }
    # If the key is (only) present in the header, drop it now
    if (length(keyHeaderIndex) > 0){
        vcf <- .dropInfoHeader(vcf = vcf, key = key)
    }
    
    # Append new header
    info(header(vcf))[key,] <- newInfoHeader
    # Append values in data
    info(vcf)[,key] <- .countGenos(
        x = geno(vcf)[["GT"]][,samples], genos = genos)

    return(vcf)
}
