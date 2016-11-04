
# vcf=ExpandedVCF ----

setMethod(
    "vepInPhenoLevel", c("ExpandedVCF"),
    function(vcf, phenoCol, level, vepCol, unique = FALSE){

        return(.vepInPhenoLevel(vcf, phenoCol, level, vepCol, unique))
    }
)

# Main method ----

# Extract VEP prediction for variants seen in a phenotype level
# vcf = ExpandedVCF
# phenoCol = character(1)
# level = character(1)
# unique = logical(1)
# facet = character(n) or NULL
.vepInPhenoLevel <- function(
    vcf, phenoCol, level, vepCol,
    unique, facet){

    # Validate input arguments
    stopifnot(is.character(phenoCol))
    stopifnot(length(phenoCol) == 1)

    stopifnot(is.character(level))
    stopifnot(length(level) == 1)

    stopifnot(is.character(vepCol))

    stopifnot(is.logical(unique))
    stopifnot(length(unique) == 1)

    stopifnot("TVTBparam" %in% names(metadata(vcf)))
    param <- metadata(vcf)[["TVTBparam"]]

    # Shortcut
    phenos <- colData(vcf)

    # Identify samples with the phenotype of interest
    samplesIdx <- which(phenos[,phenoCol] == level)

    # Identify the variants observed in phenotype level
    variantsIdx <- .variantsInSamples(vcf, samplesIdx, unique)

    # parse VEP predictions of variants in phenotype level (maybe no variant)
    csq <- parseCSQToGRanges(vcf[variantsIdx], info.key = vep(param))

    # Cannot selected columns of an empty GRanges
    return(csq[,vepCol])
}
