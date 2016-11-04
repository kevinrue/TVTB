
### phenos = list ----
## param = tSVEParam ----

setMethod(
    f = "addFrequencies",
    signature = c(vcf="ExpandedVCF", phenos="list", param="tSVEParam"),
    definition = function(
        vcf, phenos, param, ..., force = FALSE){

        param <- .override.tSVEParam(param, ...)

        .addFrequencies(
            vcf = vcf, param = param, phenos = phenos, force = force)
    }
)

### param = missing ----

setMethod(
    f = "addFrequencies",
    signature = c(vcf="ExpandedVCF", phenos="list", param="missing"),
    definition = function(
        vcf, ref, het, alt, phenos, ..., force = FALSE){

        # Use default tSVEParams
        param <- tSVEParam(genos = list(ref, het, alt))

        .addFrequencies(
            vcf = vcf, param = param, phenos = phenos, force = force)
    }
)

### phenos = character ----
## param = tSVEParam ----

setMethod(
    f = "addFrequencies",
    signature = c(vcf="ExpandedVCF", phenos="character", param="tSVEParam"),
    definition = function(
        vcf, phenos, param, ..., force = FALSE){

        param <- .override.tSVEParam(param, ...)

        # named list of all levels of given phenotypes
        phenos <- sapply(X = phenos, FUN = function(x){
            unique(colData(vcf)[,x])
        })

        .addFrequencies(
            vcf = vcf, param = param, phenos = phenos, force = force)
    }
)

### param = missing ----

setMethod(
    f = "addFrequencies",
    signature = c(vcf="ExpandedVCF", phenos="character", param="missing"),
    definition = function(
        vcf, ref, het, alt, phenos, ..., force = FALSE){

        # Use default tSVEParams
        param <- tSVEParam(genos = list(ref, het, alt))

        # named list of all levels of given phenotypes
        phenos <- sapply(X = phenos, FUN = function(x){
            unique(colData(vcf)[,x])
        })

        .addFrequencies(
            vcf = vcf, param = param, phenos = phenos, force = force)
    }
)

### phenos = missing ----
## param = tSVEParam ----

setMethod(
    f = "addFrequencies",
    signature = c(vcf="ExpandedVCF", phenos="missing", param="tSVEParam"),
    definition = function(
        vcf, param, ..., force = FALSE){

        param <- .override.tSVEParam(param, ...)

        .addFrequencies(
            vcf = vcf, param = param, phenos = list(), force = force)
    }
)

### param = missing ----

setMethod(
    f = "addFrequencies",
    signature = c(vcf="ExpandedVCF", phenos="missing", param="missing"),
    definition = function(
        vcf, ref, het, alt, ..., force = FALSE){

        # Use default tSVEParams
        param <- tSVEParam(genos = list(ref, het, alt))

        .addFrequencies(
            vcf = vcf, param = param, phenos = list(), force = force)
    }
)

# Main method ----

.checkFrequencyInfo <- function(vcf, param, phenos, force){
    # Collate all the phenotypes/levels pairs present
    phenoLevels <- unlist(lapply(
        X = 1:length(phenos),
        FUN = function(x){paste(names(phenos)[x], phenos[[x]], sep = "_")})
    )
    # Collate the INFO keys
    keySuffixes <- c(names(genos(param)), aaf(param), maf(param))
    # Deduce all INFO keys needed
    infoKeys <- as.character(unlist(sapply(
        X = phenoLevels,
        FUN = function(x){paste(x, keySuffixes, sep = "_")},
        simplify = FALSE))
    )

    matches <- match(infoKeys, colnames(info(vcf)))
    idxMatches <- matches[!is.na(matches)]

    if ((length(idxMatches) > 0))
        if (force){
            # Remove data and header
            message("Overwriting INFO fields: ", colnames(info(vcf))[matches])
            info(vcf) <- info(vcf)[,-idxMatches]
            info(header(vcf)) <- info(header(vcf))[-idxMatches,]
        } else{
            stop("INFO keys already present:", colnames(info(vcf))[matches])
        }

    return(vcf)
}

.addFrequencies <- function(
    vcf, param, phenos = list(), force = FALSE){

    if ((length(phenos) > 0) & (!force))
        vcf <- .checkFrequencyInfo(
            vcf = vcf, param = param, phenos = phenos, force = force)

    if (length(phenos) == 0){
        vcf <- .addOverallFrequencies(vcf = vcf, param = param, force = force)
    } else {
        # TODO: Parallelise by phenotype
        for (pheno in names(phenos)){
            for (level in phenos[[pheno]]){

                vcf <- .addPhenoLevelFrequencies(
                    vcf = vcf,
                    param = param,
                    pheno = pheno,
                    level = level,
                    force = force)
            }
        }
    }

    return(vcf)
}
