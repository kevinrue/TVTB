
# vcf=ExpandedVCF ----

setMethod(
    "addPhenoLevelFrequencies", c("ExpandedVCF"),
    function(vcf, pheno, level, force = FALSE){

        return(.addPhenoLevelFrequencies(vcf, pheno, level, force))
    }
)

# Checks ----

# Check type and length of inputs
.checkInputsPLF <- function(vcf, pheno, level, force){

    stopifnot("TVTBparam" %in% names(metadata(vcf)))

    stopifnot(is.character(pheno))
    stopifnot(length(pheno) == 1)

    stopifnot(is.character(level))
    stopifnot(length(level) == 1)

    stopifnot(is.logical(force))

    stopifnot(pheno %in% colnames(colData(vcf)))

    stopifnot(level %in% colData(vcf)[,pheno])

    return(TRUE)
}

# Check that none of the required INFO keys already exist, or drop them
.checkPLFInfo <- function(vcf, pheno, level, force){

    param <- metadata(vcf)[["TVTBparam"]]

    # Collate the INFO key suffixes
    keySuffixes <- suffix(param)

    # Deduce all INFO keys required
    infoKeys <- paste(pheno, level, keySuffixes, sep = "_")

    # Index of required key in INFO header and data
    matchesHeader <- na.omit(match(infoKeys, rownames(info(header(vcf)))))
    matchesData <- na.omit(match(infoKeys, colnames(info(vcf))))

    # Process data first to avoid VCF validity warning
    if (length(matchesData) > 0){
        if (force){
            # Drop data
            message(
                "Overwriting INFO keys in data:\n- ",
                paste(colnames(info(vcf))[matchesData], collapse = "\n- "))
            info(vcf) <- info(vcf)[,-matchesData, drop = FALSE]
        } else {
            # Throw an error
            stop(
                "INFO keys already present in data:\n- ",
                paste(colnames(info(vcf))[matchesData], collapse = "\n- "))
        }
    }

    # Process header last to avoid VCF validity warning
    if (length(matchesHeader) > 0){
        if (force){
            # Remove fields from data
            message(
                "Overwriting INFO keys in header:\n- ",
                paste(
                    rownames(info(header(vcf)))[matchesHeader],
                    collapse = "\n- "))
            info(header(vcf)) <- info(header(vcf))[-matchesHeader,]
        } else{
            # Remove fields from header
            stop(
                "INFO keys already present in header:\n- ",
                paste(
                    rownames(info(header(vcf)))[matchesHeader],
                    collapse = "\n- "))
        }
    }

    # Return vcf, trimmed if necessary
    return(vcf)
}

# Main method ----

# vcf=ExpandedVCF
# pheno = character(1)
# level = character(1)
# force = logical(1)
.addPhenoLevelFrequencies <- function(vcf, pheno, level, force){

    # Check type and number of relevant input arguments
    stopifnot(.checkInputsPLF(vcf, pheno, level, force))

    # Check presence of required INFO keys: drop or throw error
    vcf <- .checkPLFInfo(vcf, pheno, level, force)

    # Shortcut
    param <- metadata(vcf)[["TVTBparam"]]

    # Subset of samples associated with phenotype level
    # (avoid letting countGenos do it three times below)
    GTSubset <- geno(vcf)[["GT"]][,colData(vcf)[,pheno] == level]

    # TODO: could launch 3 parallel threads to count genotypes faster (?)
    # Count of REF, HET, ALT genotypes
    REF <- .countGenos(GTSubset, as.character(ref(param)))
    HET <- .countGenos(GTSubset, as.character(het(param)))
    ALT <- .countGenos(GTSubset, as.character(alt(param)))

    # Alternate allele frequency
    AAF <- (HET + 2 * ALT) / (2 * (REF + HET + ALT))

    # Minor allele frequency
    MAF <- mapply(
        function(ref, alt){min(ref, alt)},
        ref = 1 - AAF,
        alt = AAF) # TODO: re-introduce parallel processing

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
        # NOTE: currently,  object: 'info(VCFHeader)' must be a 3 column
        # DataFrame with names Number, Type, Description
        # VCF4.2 format mentions fields "Source" and "Version" (?)
        #Source = rep("TVTB", 5),
        #Version = rep(packageVersion("TVTB"), 5),
        Number = rep(1, 5),
        Type = c(rep("Integer", 3), rep("Float", 2)),
        Description = c(desc_REF, desc_HET, desc_ALT, desc_AAF, desc_MAF),
        row.names = paste(
            pheno,
            level,
            suffix(param)[c("ref", "het", "alt", "aaf", "maf")],
            sep = "_")
    )

    # Collate new data and column names
    newInfoData <- DataFrame(REF, HET, ALT, AAF, MAF)
    colnames(newInfoData) <- rownames(newInfoHeader)

    # Append new header fields first to avoid validity check warning
    info(header(vcf)) <- rbind(info(header(vcf)), newInfoHeader)

    # Append new data fields last to avoid validity check warning
    info(vcf) <- cbind(
        info(vcf),
        newInfoData)

    return(vcf)
}
