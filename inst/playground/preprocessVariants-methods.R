### file = TabixFile ----
## param = TVTBparam ----

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", param="TVTBparam",
        phenos="DataFrame"),
    definition = function(file, param, phenos, ...){

        param <- .override.TVTBparam(param = param, ...)

        .preprocessVariants(file = file, param = param, phenos = phenos, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", param="TVTBparam",
        phenos="data.frame"),
    definition = function(file, param, phenos, ...){

        param <- .override.TVTBparam(param = param, ...)

        .preprocessVariants(
            file = file, param = param,
            phenos = DataFrame(phenos), ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", param="TVTBparam",
        phenos="character"),
    definition = function(file, param, phenos, ...){

        param <- .override.TVTBparam(param = param, ...)

        phenosDF <- DataFrame(read.table(
            file = phenos,
            header = TRUE,
            row.names = 1,
            ...))

        .preprocessVariants(
            file = file, param = param,
            phenos = phenosDF, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", param="TVTBparam",
        phenos="missing"),
    definition = function(file, param, phenos = DataFrame(), ...){

        param <- .override.TVTBparam(param = param, ...)

        .preprocessVariants(
            file = file, param = param,
            phenos = DataFrame(), ...)
    }
)

### param = missing ----

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", param="missing",
        phenos="DataFrame"),
    definition = function(
        file,
        param = TVTBparam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos, ...){

        # Only required for vep field
        param <- TVTBparam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.TVTBparam(param = param, ...)

        .preprocessVariants(
            file = file, param = param,
            phenos = phenos, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", param="missing",
        phenos="data.frame"),
    definition = function(
        file,
        param = TVTBparam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos, ...){

        # Only required for vep field
        param <- TVTBparam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.TVTBparam(param = param, ...)

        .preprocessVariants(
            file = file, param = param,
            phenos = DataFrame(phenos), ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", param="missing",
        phenos="character"),
    definition = function(
        file,
        param = TVTBparam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos, ...){

        # Only required for vep field
        param <- TVTBparam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.TVTBparam(param = param, ...)

        phenosDF <- DataFrame(read.table(
            file = phenos,
            header = TRUE,
            row.names = 1,
            ...))

        .preprocessVariants(
            file = file, param = param,
            phenos = phenosDF, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", param="missing",
        phenos="missing"),
    definition = function(
        file,
        param = TVTBparam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos = DataFrame(), ...){

        # Only required for vep field
        param <- TVTBparam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.TVTBparam(param = param, ...)

        .preprocessVariants(
            file = file, param = param,
            phenos = DataFrame(), ...)
    }
)

### file = character ----
## param = TVTBparam ----

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", param="TVTBparam",
        phenos="DataFrame"),
    definition = function(file, param, phenos, ...){

        param <- .override.TVTBparam(param = param, ...)

        .preprocessVariants(
            file = file, param = param,
            phenos = phenos, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", param="TVTBparam",
        phenos="data.frame"),
    definition = function(file, param, phenos, ...){

        param <- .override.TVTBparam(param = param, ...)

        .preprocessVariants(
            file = file, param = param,
            phenos = DataFrame(phenos), ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", param="TVTBparam",
        phenos="character"),
    definition = function(file, param, phenos, ...){

        param <- .override.TVTBparam(param = param, ...)

        phenosDF <- DataFrame(read.table(
            file = phenos,
            header = TRUE,
            row.names = 1,
            ...))

        .preprocessVariants(
            file = file, param = param,
            phenos = phenosDF, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", param="TVTBparam",
        phenos="missing"),
    definition = function(file, param, phenos = DataFrame(), ...){

        param <- .override.TVTBparam(param = param, ...)

        .preprocessVariants(
            file = file, param = param,
            phenos = DataFrame(), ...)
    }
)

## param = missing ----

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", param="missing",
        phenos="DataFrame"),
    definition = function(
        file,
        param = TVTBparam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos, ...){

        # Only required for vep field
        param <- TVTBparam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.TVTBparam(param = param, ...)

        .preprocessVariants(
            file = file, param = param,
            phenos = phenos, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", param="missing",
        phenos="data.frame"),
    definition = function(
        file,
        param = TVTBparam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos, ...){

        # Only required for vep field
        param <- TVTBparam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.TVTBparam(param = param, ...)

        .preprocessVariants(
            file = file, param = param,
            phenos = DataFrame(phenos), ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", param="missing",
        phenos="character"),
    definition = function(
        file,
        param = TVTBparam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos, ...){

        # Only required for vep field
        param <- TVTBparam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.TVTBparam(param = param, ...)

        phenosDF <- DataFrame(read.table(
            file = phenos,
            header = TRUE,
            row.names = 1,
            ...))

        .preprocessVariants(
            file = file, param = param,
            phenos = phenosDF, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", param="missing",
        phenos="missing"),
    definition = function(
        file,
        param = TVTBparam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos = DataFrame(), ...){

        # Only required for vep field
        param <- TVTBparam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.TVTBparam(param = param, ...)

        .preprocessVariants(
            file = file, param = param,
            phenos = DataFrame(), ...)
    }
)

# Main method ----

.preprocessVariants <- function(
    file, param,
    phenos = DataFrame(), ...){

    vcfPath <- switch (class(file),
        TabixFile = path(file),
        character = file
    )

    vcfHead <- scanVcfHeader(file = file)

    if (!vep(param) %in% rownames(info(vcfHead))){
        stop(vep(param), " not found in VCF header")
    }

    svp <- ScanVcfParam(
        fixed = "ALT",
        info = vep(param),
        geno = c("GT")
    )

    if (nrow(phenos) > 0){
        if (!any(rownames(phenos) %in% samples(vcfHead))){
            stop("All samples in Phenotype file are missing from VCF header")
        }
        if (!all(rownames(phenos) %in% samples(vcfHead))){
            stop(
                "Some of the samples in Phenotype file are missing from VCF",
                " header")
        }
        # All phenoSamples in VCF: import only phenoSamples from VCF
        vcfSamples(svp) <- rownames(phenos)
    }

    if (length(ranges(param)) > 0){
        vcfWhich(svp) <- ranges(param)
    }

    message(
        "Importing VCF fields ",
        "(ALT, ",
        "FORMAT/GT ",
        "INFO/", vep(param),
        ") from ",
        length(ranges(param)), " region(s) in ",
        vcfPath,
        " ...")

    vcf <- readVcf(file = file, param = svp, ...)

    # NOTE: currently,  object: 'info(VCFHeader)' must be a 3 column DataFrame
    # with names Number, Type, Description
    # If missing, add Source and Version fields
    # if (! "Source" %in% colnames(info(header(vcf)))){
    #     info(header(vcf))[,"Source"] <- "NA"
    # }
    # if (! "Version" %in% colnames(info(header(vcf)))){
    #     info(header(vcf))[,"Version"] <- "NA"
    # }

    # Attach phenotypes information if given
    if (nrow(phenos) > 0){
        colData(vcf) <- phenos
    }

    # Expand to bi-allelic records
    vcf <- expand(x = vcf, row.names = TRUE)

    # If any multi-allelic record is present in the file
    #if (any(lengths(mcols(vcf)[,"ALT"]) > 1)){
        # Disambiguate their row.names using, appending "_ALT"
    rownames(vcf) <- paste(rownames(vcf), mcols(vcf)[,"ALT"], sep = "_")
    #}

    # Header of all INFO fields were imported but not all data fields,
    # remove headers of fields not imported
    info(header(vcf)) <- subset(
        x = info(header(vcf)),
        subset = rownames(info(header(vcf))) %in% colnames(info(vcf)))

    return(vcf)
}
