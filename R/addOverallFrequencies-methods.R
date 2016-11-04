
# vcf=ExpandedVCF ----

setMethod(
    "addOverallFrequencies", c("ExpandedVCF"),
    function(vcf, force = FALSE){

        return(.addOverallFrequencies(vcf, force))
    }
)

# Main method ----

# vcf = ExpandedVCF
# force = logical(1)
.checkOverallInfo <- function(vcf, force){

    stopifnot("TVTBparam" %in% names(metadata(vcf)))
    param <- metadata(vcf)[["TVTBparam"]]

    infoKeys <- c(names(genos(param)), aaf(param), maf(param))

    matchesHeader <- na.omit(match(infoKeys, rownames(info(header(vcf)))))
    matchesData <- na.omit(match(infoKeys, colnames(info(vcf))))

    # Process data first to avoid VCF validity warning
    if ((length(matchesData) > 0)){
        if (force){
            # Remove fields from data
            message(
                "Overwriting INFO keys in data:\n- ",
                paste(colnames(info(vcf))[matchesData], collapse = "\n- "))
            info(vcf) <- info(vcf)[,-matchesData, drop = FALSE]
        } else{
            stop(
                "INFO keys already present in data: ",
                paste(colnames(info(vcf))[matchesData], collapse = "\n- "))
        }
    }

    # Process header last to avoid VCF validity warning
    if ((length(matchesHeader) > 0)){
        if (force){
            # Remove fields from header
            message(
                "Overwriting INFO keys in header:\n- ",
                paste(
                    rownames(info(header(vcf)))[matchesHeader],
                    collapse = "\n- "))
            info(header(vcf)) <- info(header(vcf))[-matchesHeader,]
        } else{
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

# vcf = ExpandedVCF
# force = logical(1)
.addOverallFrequencies <- function(vcf, force){

    # Check presence of required INFO keys: drop or throw error
    vcf <- .checkOverallInfo(vcf, force)

    # Shortcut
    GT <- geno(vcf)[["GT"]]
    param <- metadata(vcf)[["TVTBparam"]]

    # TODO: could launch 3 parallel threads to count genotypes
    # Count of REF, HET, ALT genotypes
    REF <- .countGenos(GT, ref(param))
    HET <- .countGenos(GT, het(param))
    ALT <- .countGenos(GT, alt(param))

    # Alternate allele frequency
    AAF <- (HET + 2 * ALT) / (2 * (REF + HET + ALT))

    # Minor allele frequency
    MAF <- bpmapply(
        function(ref, alt){min(ref, alt)},
        ref = 1 - AAF,
        alt = AAF,
        BPPARAM = bp(param))

    # Collate new headers
    newInfoHeader <- DataFrame(
        Number = rep(1, 5),
        Type = c(rep("Integer", 3), rep("Float", 2)),
        Description = c(
            "Number of homozygote reference genotypes",
            "Number of heterozygous genotypes",
            "Number of homozygote alternate genotypes",
            "Alternate allele frequency",
            "Minor allele frequency"
        ),
        # NOTE: currently,  object: 'info(VCFHeader)' must be a 3 column
        # DataFrame with names Number, Type, Description
        #Source = rep("TVTB", 5),
        #Version = rep(packageVersion("TVTB"), 5),
        row.names = suffix(param)[c("ref", "het", "alt", "aaf", "maf")]
    )

    # Collate new data
    newInfoData <- DataFrame(REF, HET, ALT, AAF, MAF)
    colnames(newInfoData) <- rownames(newInfoHeader)

    # Append new header fields
    info(header(vcf)) <- rbind(info(header(vcf)), newInfoHeader)

    # Append new data fields
    info(vcf) <- cbind(
        info(vcf),
        newInfoData)

    return(vcf)
}
