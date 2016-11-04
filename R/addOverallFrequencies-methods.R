
## param = TVTBparam ----

setMethod(
    f = "addOverallFrequencies",
    signature = c(vcf="ExpandedVCF", param="TVTBparam"),
    definition = function(
        vcf, param, ..., force = FALSE){

        param <- .override.TVTBparam(param = param, ...)

        .addOverallFrequencies(
            vcf = vcf, param = param, force = force)
    }
)

# Main method ----

.checkOverallInfo <- function(vcf, param, force){

    infoKeys <- c(names(genos(param)), aaf(param), maf(param))

    matchesHeader <- na.omit(match(infoKeys, rownames(info(header(vcf)))))
    matchesData <- na.omit(match(infoKeys, colnames(info(vcf))))

    # Data first to avoid validity warning
    if ((length(matchesData) > 0))
        if (force){
            # Remove data and header
            message(
                "Overwriting INFO keys in data:\n- ",
                paste(colnames(info(vcf))[matchesData], collapse = "\n- "))
            info(vcf) <- info(vcf)[,-matchesData, drop = FALSE]
        } else{
            stop(
                "INFO keys already present in data: ",
                paste(colnames(info(vcf))[matchesData], collapse = "\n- "))
        }

    if ((length(matchesHeader) > 0))
        if (force){
            # Remove data and header
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

    return(vcf)

}

.addOverallFrequencies <- function(vcf, param, force = FALSE){

    vcf <- .checkOverallInfo(vcf = vcf, param = param, force = force)

    GT <- geno(vcf)[["GT"]]
    # TODO: could launch 3 parallel threads to count genotypes
    # Count of REF, HET, ALT genotypes
    REF <- .countGenos(x = GT, genos = hRef(param))
    HET <- .countGenos(x = GT, genos = het(param))
    ALT <- .countGenos(x = GT, genos = hAlt(param))

    # Alternate allele frequency
    AAF <- (HET + 2 * ALT) / (2 * (REF + HET + ALT))

    # Minor allele frequency
    MAF <- bpmapply(
        FUN = function(ref, alt){min(ref, alt)},
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
        row.names = c(names(genos(param)), aaf(param), maf(param))
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
