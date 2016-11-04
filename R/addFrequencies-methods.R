
### phenos = list ----
## param = TVTBparam ----

setMethod(
    f = "addFrequencies",
    signature = c(vcf="ExpandedVCF", phenos="list", param="TVTBparam"),
    definition = function(
        vcf, phenos, param, ..., force = FALSE){

        param <- .override.TVTBparam(param, ...)

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

        # Use default TVTBparam
        param <- TVTBparam(genos = list(ref, het, alt))

        .addFrequencies(
            vcf = vcf, param = param, phenos = phenos, force = force)
    }
)

### phenos = character ----
## param = TVTBparam ----

setMethod(
    f = "addFrequencies",
    signature = c(vcf="ExpandedVCF", phenos="character", param="TVTBparam"),
    definition = function(
        vcf, phenos, param, ..., force = FALSE){

        param <- .override.TVTBparam(param, ...)

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

        # Use default TVTBparam
        param <- TVTBparam(genos = list(ref, het, alt))

        # named list of all levels of given phenotypes
        phenos <- sapply(X = phenos, FUN = function(x){
            unique(colData(vcf)[,x])
        })

        .addFrequencies(
            vcf = vcf, param = param, phenos = phenos, force = force)
    }
)

### phenos = missing ----
## param = TVTBparam ----

setMethod(
    f = "addFrequencies",
    signature = c(vcf="ExpandedVCF", phenos="missing", param="TVTBparam"),
    definition = function(
        vcf, param, ..., force = FALSE){

        param <- .override.TVTBparam(param, ...)

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

        # Use default TVTBparam
        param <- TVTBparam(genos = list(ref, het, alt))

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

    matches <- as.numeric(na.omit(match(infoKeys, colnames(info(vcf)))))

    if ((length(matches) > 0))
        stop(
            "INFO keys already present: ",
            paste(colnames(info(vcf))[matches], sep = ", "))

    return(TRUE)
}

.addFrequencies <- function(
    vcf, param, phenos = list(), force = FALSE){

    if ((length(phenos) > 0) & (!force))
        .checkFrequencyInfo(
            vcf = vcf, param = param, phenos = phenos)

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
