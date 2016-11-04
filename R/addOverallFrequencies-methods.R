
## param = tSVEParam ----

setMethod(
    f = "addOverallFrequencies",
    signature = c(vcf="ExpandedVCF", param="tSVEParam"),
    definition = function(
        vcf, param, ..., force = FALSE){

        param <- .override.tSVEParam(param = param, ...)

        .addOverallFrequencies(
            vcf = vcf, param = param, force = force)
    }
)

## param = missing ----

setMethod(
    f = "addOverallFrequencies",
    signature = c(vcf="ExpandedVCF", param="missing"),
    definition = function(
        vcf, ref, het, alt, ..., force = FALSE){

        # Use default tSVEParams
        param <- tSVEParam(genos = list(ref, het, alt))

        .addOverallFrequencies(
            vcf = vcf, param = param, force = force)
    }
)

# Main method ----

.checkOverallInfo <- function(vcf, param, force){

    infoKeys <- c(names(genos(param)), aaf(param), maf(param))

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

.addOverallFrequencies <- function(vcf, param, force = FALSE){

    .checkOverallInfo(vcf = vcf, param = param, force = force)

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
        #Source = rep("tSVE", 5),
        #Version = rep(packageVersion("tSVE"), 5),
        row.names = c(names(genos(param)), aaf(param), maf(param))
    )

    # Collate new data
    newInfoData <- DataFrame(
        REF = REF,
        HET = HET,
        ALT = ALT,
        AAF = AAF,
        MAF = MAF
    )

    # Append new header fields
    info(header(vcf)) <- rbind(info(header(vcf)), newInfoHeader)

    # Append new data fields
    info(vcf) <- cbind(
        info(vcf),
        newInfoData)

    return(vcf)
}
