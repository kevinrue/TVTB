### param = TVTBparam ----
## phenotypes = DataFrame ----

setMethod(
    f = "vepInPhenoLevel",
    signature = c(vcf="ExpandedVCF", param="TVTBparam"),
    definition = function(
        vcf, phenoCol, level, vepCol, param = NULL, ...,
        unique = FALSE, facet = NULL){

        param <- .override.TVTBparam(param, ...)

        .vepInPhenoLevel(
            vcf = vcf, phenoCol = phenoCol, level = level, vepCol = vepCol,
            param = param,
            unique = unique, facet = facet)}
)

### param = missing ----
## phenotypes = DataFrame ----

setMethod(
    f = "vepInPhenoLevel",
    signature = c(vcf="ExpandedVCF", param="missing"),
    definition = function(
        vcf, phenoCol, level, vepCol, alts, param = NULL, ...,
        unique = FALSE, facet = NULL){

        if (length(alts) < 2)
            stop(
                "length(alts) must be >= 2: ",
                "Heterozygote and Homozygote alternate genotypes")
        # ref will not be used
        param <- TVTBparam(genos = list("", alts[1], alts[2:length(alts)]))

        param <- .override.TVTBparam(param, ...)

        .vepInPhenoLevel(
            vcf = vcf, phenoCol = phenoCol, level = level, vepCol = vepCol,
            param = param,
            unique = unique, facet = facet)}
)

.vepInPhenoLevel <- function(
    vcf, phenoCol, level, vepCol, param,
    unique = FALSE, facet = NULL){

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
