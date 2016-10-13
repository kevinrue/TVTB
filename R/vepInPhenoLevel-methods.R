
# vcf=ExpandedVCF ----

setMethod(
    "vepInPhenoLevel", c("ExpandedVCF"),
    function(
        vcf, phenoCol, level, vepCol,
        unique = FALSE, facet = NULL){

        return(.vepInPhenoLevel(
            vcf, phenoCol, level, vepCol,
            unique, facet))
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
    stopifnot(length(vepCol) == 1)

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

    # Extract requested predictions
    if (length(csq) > 0){
        # Keep the desired consequence(+facet) for those variants
        vepFacet <- as.data.frame(mcols(csq[,c(vepCol, facet)]))
        # Set column names (facet may be NULL)
        colnames(vepFacet) <- c(vepCol, facet)
        # Return the formatted table
        return(vepFacet)
    }

    # If no VEP predictions/variants in phenotype level: empty data.frame
    return(data.frame())
}
