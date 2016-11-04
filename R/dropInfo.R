#vcf, key = NULL, slot = "header"|"data"|"both"

setMethod(
    f = "dropInfo",
    signature = c(vcf="VCF"),
    definition = function(vcf, key = NULL, slot = "both"){
        .dropInfo(vcf, key = key, slot = slot)
    }
)

# main method ----

.dropInfo <- function(vcf, key = NULL, slot = "both"){
    switch(
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
}

.dropInfoHeader <- function(vcf, key = NULL){
    if (!is.null(key)){
        matches <- match(key, rownames(info(header(vcf))))
        missingIdx <- which(is.na(matches))
        if (length(missingIdx) > 0)
            message(
                "key not in header: ",
                paste(key[missingIdx], collapse = ", "))
        idx <- na.omit(matches)
    } else {
        # remove header without data
        idx <- which(!rownames(info(header(vcf))) %in% colnames(info(vcf)))
    }

    if (length(idx) > 0)
        info(header(vcf)) <- info(header(vcf))[-idx,]

    return(vcf)
}

.dropInfoData <- function(vcf, key = NULL){
    if (!is.null(key)){
        matches <- match(key, colnames(info(vcf)))
        missingIdx <- which(is.na(matches))
        if (length(missingIdx) > 0)
            message(
                "key not in data: ",
                paste(key[missingIdx], collapse = ", "))
        idx <- na.omit(matches)
    } else {
        # remove header without data
        idx <- which(!colnames(info(vcf)) %in% rownames(info(header(vcf))))
    }

    if (length(idx) > 0)
        info(vcf) <- info(vcf)[,-idx, drop = FALSE]

    return(vcf)
}
