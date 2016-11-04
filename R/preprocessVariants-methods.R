#### regions = GRanges ----
### file = TabixFile ----
## param = tSVEParam ----

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", regions="GRanges", param="tSVEParam",
        phenos="DataFrame"),
    definition = function(file, regions, param, phenos, ...){

        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = phenos, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", regions="GRanges", param="tSVEParam",
        phenos="data.frame"),
    definition = function(file, regions, param, phenos, ...){

        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = DataFrame(phenos), ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", regions="GRanges", param="tSVEParam",
        phenos="character"),
    definition = function(file, regions, param, phenos, ...){

        param <- .override.tSVEParam(param = param, ...)

        phenosDF <- DataFrame(read.table(
            file = phenos,
            header = TRUE,
            row.names = 1,
            ...))

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = phenosDF, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", regions="GRanges", param="tSVEParam",
        phenos="missing"),
    definition = function(file, regions, param, phenos = DataFrame(), ...){

        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = DataFrame(), ...)
    }
)

### param = missing ----

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", regions="GRanges", param="missing",
        phenos="DataFrame"),
    definition = function(
        file, regions,
        param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos, ...){

        # Only required for vep field
        param <- tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = phenos, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", regions="GRanges", param="missing",
        phenos="data.frame"),
    definition = function(
        file, regions,
        param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos, ...){

        # Only required for vep field
        param <- tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = DataFrame(phenos), ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", regions="GRanges", param="missing",
        phenos="character"),
    definition = function(
        file, regions,
        param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos, ...){

        # Only required for vep field
        param <- tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.tSVEParam(param = param, ...)

        phenosDF <- DataFrame(read.table(
            file = phenos,
            header = TRUE,
            row.names = 1,
            ...))

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = phenosDF, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", regions="GRanges", param="missing",
        phenos="missing"),
    definition = function(
        file, regions,
        param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos = DataFrame(), ...){

        # Only required for vep field
        param <- tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = DataFrame(), ...)
    }
)

### file = character ----
## param = tSVEParam ----

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", regions="GRanges", param="tSVEParam",
        phenos="DataFrame"),
    definition = function(file, regions, param, phenos, ...){

        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = phenos, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", regions="GRanges", param="tSVEParam",
        phenos="data.frame"),
    definition = function(file, regions, param, phenos, ...){

        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = DataFrame(phenos), ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", regions="GRanges", param="tSVEParam",
        phenos="character"),
    definition = function(file, regions, param, phenos, ...){

        param <- .override.tSVEParam(param = param, ...)

        phenosDF <- DataFrame(read.table(
            file = phenos,
            header = TRUE,
            row.names = 1,
            ...))

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = phenosDF, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", regions="GRanges", param="tSVEParam",
        phenos="missing"),
    definition = function(file, regions, param, phenos = DataFrame(), ...){

        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = DataFrame(), ...)
    }
)

## param = missing ----

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", regions="GRanges", param="missing",
        phenos="DataFrame"),
    definition = function(
        file, regions,
        param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos, ...){

        # Only required for vep field
        param <- tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = phenos, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", regions="GRanges", param="missing",
        phenos="data.frame"),
    definition = function(
        file, regions,
        param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos, ...){

        # Only required for vep field
        param <- tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = DataFrame(phenos), ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", regions="GRanges", param="missing",
        phenos="character"),
    definition = function(
        file, regions,
        param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos, ...){

        # Only required for vep field
        param <- tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.tSVEParam(param = param, ...)

        phenosDF <- DataFrame(read.table(
            file = phenos,
            header = TRUE,
            row.names = 1,
            ...))

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = phenosDF, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", regions="GRanges", param="missing",
        phenos="missing"),
    definition = function(
        file, regions,
        param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos = DataFrame(), ...){

        # Only required for vep field
        param <- tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = DataFrame(), ...)
    }
)

#### regions = missing ----
### file = TabixFile ----
## param = tSVEParam ----

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", regions="missing", param="tSVEParam",
        phenos="DataFrame"),
    definition = function(file, regions, param, phenos, ...){

        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = phenos, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", regions="missing", param="tSVEParam",
        phenos="data.frame"),
    definition = function(file, regions, param, phenos, ...){

        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = DataFrame(phenos), ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", regions="missing", param="tSVEParam",
        phenos="character"),
    definition = function(file, regions, param, phenos, ...){

        param <- .override.tSVEParam(param = param, ...)

        phenosDF <- DataFrame(read.table(
            file = phenos,
            header = TRUE,
            row.names = 1,
            ...))

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = phenosDF, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", regions="missing", param="tSVEParam",
        phenos="missing"),
    definition = function(file, regions, param, phenos = DataFrame(), ...){

        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = DataFrame(), ...)
    }
)

### param = missing ----

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", regions="missing", param="missing",
        phenos="DataFrame"),
    definition = function(
        file, regions,
        param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos, ...){

        # Only required for vep field
        param <- tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = phenos, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", regions="missing", param="missing",
        phenos="data.frame"),
    definition = function(
        file, regions,
        param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos, ...){

        # Only required for vep field
        param <- tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = DataFrame(phenos), ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", regions="missing", param="missing",
        phenos="character"),
    definition = function(
        file, regions,
        param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos, ...){

        # Only required for vep field
        param <- tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.tSVEParam(param = param, ...)

        phenosDF <- DataFrame(read.table(
            file = phenos,
            header = TRUE,
            row.names = 1,
            ...))

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = phenosDF, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="TabixFile", regions="missing", param="missing",
        phenos="missing"),
    definition = function(
        file, regions,
        param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos = DataFrame(), ...){

        # Only required for vep field
        param <- tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = DataFrame(), ...)
    }
)

### file = character ----
## param = tSVEParam ----

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", regions="missing", param="tSVEParam",
        phenos="DataFrame"),
    definition = function(file, regions, param, phenos, ...){

        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = phenos, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", regions="missing", param="tSVEParam",
        phenos="data.frame"),
    definition = function(file, regions, param, phenos, ...){

        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = DataFrame(phenos), ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", regions="missing", param="tSVEParam",
        phenos="character"),
    definition = function(file, regions, param, phenos, ...){

        param <- .override.tSVEParam(param = param, ...)

        phenosDF <- DataFrame(read.table(
            file = phenos,
            header = TRUE,
            row.names = 1,
            ...))

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = phenosDF, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", regions="missing", param="tSVEParam",
        phenos="missing"),
    definition = function(file, regions, param, phenos = DataFrame(), ...){

        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = DataFrame(), ...)
    }
)

## param = missing ----

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", regions="missing", param="missing",
        phenos="DataFrame"),
    definition = function(
        file, regions,
        param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos, ...){

        # Only required for vep field
        param <- tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = phenos, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", regions="missing", param="missing",
        phenos="data.frame"),
    definition = function(
        file, regions,
        param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos, ...){

        # Only required for vep field
        param <- tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = DataFrame(phenos), ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", regions="missing", param="missing",
        phenos="character"),
    definition = function(
        file, regions,
        param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos, ...){

        # Only required for vep field
        param <- tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.tSVEParam(param = param, ...)

        phenosDF <- DataFrame(read.table(
            file = phenos,
            header = TRUE,
            row.names = 1,
            ...))

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = phenosDF, ...)
    }
)

setMethod(
    f = "preprocessVariants",
    signature = c(
        file="character", regions="missing", param="missing",
        phenos="missing"),
    definition = function(
        file, regions,
        param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
        phenos = DataFrame(), ...){

        # Only required for vep field
        param <- tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.tSVEParam(param = param, ...)

        .preprocessVariants(
            file = file, regions = regions, param = param,
            phenos = DataFrame(), ...)
    }
)

# Main method ----

.preprocessVariants <- function(
    file, regions, param,
    phenos = DataFrame(), ...){

    vcfPath <- switch (class(file),
        TabixFile = path(file),
        character = file
    )

    vcfHead <- scanVcfHeader(file = file)

    if (!vep(param) %in% rownames(info(vcfHead))){
        stop(vep(param), " not found in VCF header")
    }

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
        phenoSamples <- rownames(phenos)
    } else {
        # No phenotypes: import all samples
        phenoSamples <- character()
    }

    svp <- ScanVcfParam(
        fixed = "ALT",
        info = vep(param),
        geno = c("GT"),
        samples = phenoSamples,
        which = regions
    )

    message(
        "Importing VCF fields ",
        "(ALT, ",
        "FORMAT/GT ",
        "INFO/", vep(param),
        ") from ",
        length(regions), " region(s) in ",
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
    if (any(lengths(mcols(vcf)[,"ALT"]) > 1)){
        # Disambiguate their row.names using, appending "_ALT"
        rownames(vcf) <- paste(rownames(vcf), mcols(vcf)[,"ALT"], sep = "_")
    }

    # Header of all INFO fields were imported but not all data fields,
    # remove headers of fields not imported
    extraInfo <- which(!rownames(info(header(vcf))) %in% colnames(info(vcf)))
    info(header(vcf)) <- info(header(vcf))[-extraInfo,]

    return(vcf)
}
