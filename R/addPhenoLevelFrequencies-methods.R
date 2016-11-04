
## param = TVTBparam ----

setMethod(
    f = "addPhenoLevelFrequencies",
    signature = c(vcf="ExpandedVCF", param="TVTBparam"),
    definition = function(
        vcf, pheno, level, param, ..., force = FALSE){

        param <- .override.TVTBparam(param = param, ...)

        .addPhenoLevelFrequencies(
            vcf = vcf, param = param, pheno = pheno, level = level)
    }
)

## param = missing ----

setMethod(
    f = "addPhenoLevelFrequencies",
    signature = c(vcf="ExpandedVCF", param="missing"),
    definition = function(
        vcf, pheno, level, ref, het, alt, ..., force = FALSE){

        # Use default TVTBparam
        param <- TVTBparam(genos = list(ref, het, alt))

        .addPhenoLevelFrequencies(
            vcf = vcf, param = param, pheno = pheno, level = level)
    }
)

# Main method ----

.checkInputsPLF <- function(vcf, param, pheno, level, force){

    stopifnot(length(pheno) == 1)
    stopifnot(length(level) == 1)
    stopifnot(is.logical(force))

    if (!pheno %in% colnames(colData(vcf))){
        stop(pheno, " is not a valid colnames in colData(vcf)")
    }

    if (!level %in% colData(vcf)[,pheno]){
        stop(pheno, " is not a valid value in colData(vcf)[,pheno]")
    }

    return(TRUE)
}

.checkPLFInfo <- function(vcf, param, pheno, level, force){

    # Collate the INFO keys
    keySuffixes <- c(names(genos(param)), aaf(param), maf(param))

    # Deduce all INFO keys needed
    infoKeys <- paste(pheno, level, keySuffixes, sep = "_")

    matches <- as.numeric(na.omit(match(infoKeys, colnames(info(vcf)))))

    if ((length(matches) > 0))
        if (force){
            # Remove data and header
            message("Overwriting INFO fields: ", colnames(info(vcf))[matches])
            info(vcf) <- info(vcf)[,-matches]
            info(header(vcf)) <- info(header(vcf))[-matches,]
        } else{
            stop(
                "INFO keys already present:",
                paste(colnames(info(vcf))[matches], sep = ", "))
        }

    return(vcf)
}

.addPhenoLevelFrequencies <- function(vcf, param, pheno, level, force = FALSE){

    .checkInputsPLF(
        vcf = vcf, param = param, pheno = pheno, level = level, force = force)

    vcf <- .checkPLFInfo(
        vcf = vcf, param = param, pheno = pheno, level = level, force = force)

    # Subset of samples associated with phenotype level
    # (avoid letting countGenos do it three times below)
    GTSubset <- geno(vcf)[["GT"]][,colData(vcf)[,pheno] == level]

    # TODO: could launch 3 parallel threads to count genotypes
    # Count of REF, HET, ALT genotypes
    REF <- .countGenos(x = GTSubset, genos = as.character(hRef(param)))
    HET <- .countGenos(x = GTSubset, genos = as.character(het(param)))
    ALT <- .countGenos(x = GTSubset, genos = as.character(hAlt(param)))

    # Alternate allele frequency
    AAF <- (HET + 2 * ALT) / (2 * (REF + HET + ALT))

    # Minor allele frequency
    MAF <- bpmapply(
        FUN = function(ref, alt){min(ref, alt)},
        ref = 1 - AAF,
        alt = AAF,
        BPPARAM = bp(param))

    # Collate description of new headers
    decription_suffix <- paste0(
        "in phenotype \"", pheno, "\", level \"", level, "\"")
    desc_REF <- paste(
        "Count of homozygote reference genotypes", decription_suffix)
    desc_HET <- paste(
        "Count of heterozygous genotypes", decription_suffix)
    desc_ALT <- paste(
        "Count of homozygote alternate genotypes", decription_suffix)
    desc_AAF <- paste(
        "Alternate allele frequency", decription_suffix)
    desc_MAF <- paste(
        "Minor allele frequency", decription_suffix)

    # Collate new headers
    newInfoHeader <- DataFrame(
        Number = rep(1, 5),
        Type = c(rep("Integer", 3), rep("Float", 2)),
        Description = c(desc_REF, desc_HET, desc_ALT, desc_AAF, desc_MAF),
        # NOTE: currently,  object: 'info(VCFHeader)' must be a 3 column
        # DataFrame with names Number, Type, Description
        #Source = rep("TVTB", 5),
        #Version = rep(packageVersion("TVTB"), 5),
        row.names = paste(
            pheno,
            level,
            c(names(genos(param)), aaf(param), maf(param)),
            sep = "_")
    )

    # Collate new data
    newInfoData <- DataFrame(
        REF, HET, ALT, AAF, MAF
    )
    colnames(newInfoData) <- rownames(newInfoHeader)

    # Append new header fields
    info(header(vcf)) <- rbind(info(header(vcf)), newInfoHeader)

    # Append new data fields
    info(vcf) <- cbind(
        info(vcf),
        newInfoData)

    return(vcf)
}
