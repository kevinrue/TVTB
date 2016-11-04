
# vcf=VCF ----

setMethod(
    "dropInfo", c("VCF"),
    function(vcf, key = NULL, slot = "both"){
        return(.dropInfo(vcf, key, slot))
    }
)

# Header method ----

# Drop field(s) from the INFO header
# vcf = VCF
# key = character(n) or NULL
.dropInfoHeader <- function(vcf, key){

    # Key(s) supplied: identify those present in header
    if (!is.null(key)){
        matches <- match(key, rownames(info(header(vcf))))
        missingIdx <- which(is.na(matches))
        if (length(missingIdx) > 0){
            message(
                "key(s) not in header: ",
                paste(key[missingIdx], collapse = ", "))
        }
        idx <- na.omit(matches)
    } else {
        # No key supplied: identify all header keys absent from data
        idx <- which(!rownames(info(header(vcf))) %in% colnames(info(vcf)))
    }

    # Drop relevant keys from header, if any
    if (length(idx) > 0){
        info(header(vcf)) <- info(header(vcf))[-idx,]
    }

    return(vcf)
}

# Data method ----

# Drop field(s) from the INFO data
# vcf = VCF
# key = character(n) or NULL
.dropInfoData <- function(vcf, key){

    # Key(s) supplied: identify those present in data
    if (!is.null(key)){
        matches <- match(key, colnames(info(vcf)))
        missingIdx <- which(is.na(matches))
        if (length(missingIdx) > 0)
            message(
                "key not in data: ",
                paste(key[missingIdx], collapse = ", "))
        idx <- na.omit(matches)
    } else {
        # No key supplied: identify all data keys absent from header
        idx <- which(!colnames(info(vcf)) %in% rownames(info(header(vcf))))
    }

    # Drop relevant keys from data, if any
    if (length(idx) > 0){
        info(vcf) <- info(vcf)[,-idx, drop = FALSE]
    }

    return(vcf)
}

# Main method ----

# Drop INFO fields from header and/or data
# vcf = VCF
# key = character(n) or NULL
# slot = "header", "data", or "both"
.dropInfo <- function(vcf, key, slot){

    vcf <- switch(
        slot,
        header = .dropInfoHeader(vcf, key),
        data = .dropInfoData(vcf, key),
        both = {
            # drop data first, to avoid validity warning
            vcf <- .dropInfoData(vcf, key)
            .dropInfoHeader(vcf, key)
        },
        stop("invalid slot argument")
    )

    return(vcf)
}


