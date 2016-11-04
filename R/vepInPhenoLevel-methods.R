
## phenotypes = DataFrame ----

setMethod(
    f = "vepInPhenoLevel",
    signature = c(vcf="ExpandedVCF", param="TVTBparam"),
    definition = function(
        vcf, phenoCol, level, vepCol, param, ...,
        unique = FALSE, facet = NULL){

        param <- .override.TVTBparam(param, ...)

        .vepInPhenoLevel(
            vcf = vcf, phenoCol = phenoCol, level = level, vepCol = vepCol,
            param = param,
            unique = unique, facet = facet)}
)


.vepInPhenoLevel <- function(
    vcf, phenoCol, level, vepCol, param,
    unique = FALSE, facet = NULL){

    stopifnot(identical(x = length(phenoCol), y = as.integer(1)))
    stopifnot(identical(x = length(level), y = as.integer(1)))
    
    phenos <- colData(vcf)

    # Identify samples with the phenotype of interest
    samplesIdx <- which(phenos[,phenoCol] == level)

    # Identify the variants observed in phenotype
    variantsIdx <- .variantsInSamples(
        vcf = vcf,
        samples = samplesIdx,
        param = param,
        unique = unique)

    # parse VEP predictions of those variants
    csq <- parseCSQToGRanges(
        x = vcf[variantsIdx], info.key = vep(param))

    if (length(csq) == 0)
        return(data.frame())
    # Keep the desired consequence(+facet) for those variants
    vepFacet <- as.data.frame(mcols(csq[,c(vepCol, facet)]))

    colnames(vepFacet) <- c(vepCol, facet)

    return(vepFacet)
}
