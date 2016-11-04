
# vcf=ExpandedVCF ----

setMethod(
    "addCountGenos", c("ExpandedVCF"),
    function(
        vcf, genos, key, description, samples = 1:ncol(vcf), force = FALSE){

        return(.addCountGenos(vcf, genos, key, description, samples, force))
    }
)

# Main method ----

# vcf=ExpandedVCF
# genos = character(n)
# key = character(1)
# description = character(1)
# samples = <indices or names that identify samples>
# force = logical(1)
.addCountGenos <- function(
    vcf, genos, key, description, samples, force){

    # Validate relevant inputs
    stopifnot(is.logical(force))

    stopifnot(is.character(key))
    stopifnot(length(key) == 1)

    stopifnot(is.character(description))
    stopifnot(length(description) == 1)

    # Prepare the new INFO header
    newInfoHeader <- DataFrame(
        Number = 1, Type = "Integer",
        Description = description,
        row.names = key)

    # Index of required key in INFO header and data
    keyDataIndex <- which(colnames(info(vcf)) == key)
    keyHeaderIndex <- which(rownames(info(header(vcf))) == key)

    # Check: warn/stop if key present in INFO data
    if (length(keyDataIndex) > 0){
        # If user wishes to overwrite the field
        if (force){
            # Drop the data first to avoid validity warning
            message(key, " overwritten.")
            vcf <- .dropInfoData(vcf, key)
        } else {
            # If user didn't ask to over-write the key, throw an error
            stop(key, " key already present in INFO.")
        }
    }

    # If the key is present in the header, drop it now (avoid warning)
    # NOTE: force is not required if the key is *only* present in the header
    if (length(keyHeaderIndex) > 0){
        vcf <- .dropInfoHeader(vcf, key)
    }

    # Append new header first (avoid validity warning)
    info(header(vcf))[key,] <- newInfoHeader
    # Append values in data
    info(vcf)[,key] <- .countGenos(geno(vcf)[["GT"]][,samples], genos)

    # Return updated ExpandedVCF object
    return(vcf)
}
