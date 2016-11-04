#### ByPhenotype() ----

## param = TVTBparam
setMethod(
    f = "densityVepByPhenotype",
    signature = c(vcf="ExpandedVCF", param="TVTBparam"),
    definition = function(
        vcf, phenoCol, vepCol, param, ..., filter = VcfFilterRules(),
        unique = FALSE, facet = NULL, plot = FALSE, popFreq = FALSE){

        param <- .override.TVTBparam(param, ...)

        .densityVepByPhenotype(
            vcf = vcf, phenoCol = phenoCol, vepCol = vepCol, param = param,
            filter = filter,
            unique = unique, facet = facet, plot = plot, popFreq = popFreq)
    }
)

## param = missing
setMethod(
    f = "densityVepByPhenotype",
    signature = c(vcf="ExpandedVCF", param="missing"),
    definition = function(
        vcf, phenoCol, vepCol, alts, param = NULL, ...,
        filter = VcfFilterRules(),
        unique = FALSE, facet = NULL, plot = FALSE, popFreq = FALSE){

        .checkAlts(alts)

        # ref will not be used
        param <- TVTBparam(genos = list("", alts[1], alts[2:length(alts)]))
        # Additional parameters overriden
        param <- .override.TVTBparam(param, ...)

        .densityVepByPhenotype(
            vcf = vcf, phenoCol = phenoCol, vepCol = vepCol, param = param,
            filter = filter,
            unique = unique, facet = facet, plot = plot, popFreq = popFreq)
    }
)

#### InPhenoLevel() ----

## param = TVTBparam
setMethod(
    f = "densityVepInPhenoLevel",
    signature = c(vcf="ExpandedVCF", param="TVTBparam"),
    definition = function(
        level, vcf, phenoCol, vepCol, param, ..., filter = VcfFilterRules(),
        unique = FALSE, facet = NULL, plot = FALSE, popFreq = FALSE){

        param <- .override.TVTBparam(param, ...)

        .densityVepInPhenoLevel(
            level = level, vcf = vcf, phenoCol = phenoCol, vepCol = vepCol,
            param = param, filter = filter,
            unique = unique, facet = facet, plot = plot, popFreq = popFreq)
    }
)

## param = missing
setMethod(
    f = "densityVepInPhenoLevel",
    signature = c(vcf="ExpandedVCF", param="missing"),
    definition = function(
        level, vcf, phenoCol, vepCol, alts, param = NULL, ...,
        filter = VcfFilterRules(),
        unique = FALSE, facet = NULL, plot = FALSE, popFreq = FALSE){

        .checkAlts(alts)

        # ref will not be used
        param <- TVTBparam(genos = list("", alts[1], alts[2:length(alts)]))
        # Additional parameters overriden
        param <- .override.TVTBparam(param, ...)

        .densityVepInPhenoLevel(
            level = level, vcf = vcf, phenoCol = phenoCol, vepCol = vepCol,
            param = param, filter = filter,
            unique = unique, facet = facet, plot = plot, popFreq = popFreq)
    }
)

# Helpers ----

.checkAlts <- function(alts){
    if (length(alts) < 2)
        stop(
            "length(alts) must be >= 2: ",
            "Heterozygote and Homozygote alternate genotypes")

    return(TRUE)
}

# Main methods ----

.densityVepByPhenotype <- function(
    vcf, phenoCol, vepCol, param, filter = VcfFilterRules(),
    unique = FALSE, facet = NULL, plot = FALSE, popFreq = FALSE){

    phenos <- colData(vcf)
    pLevels <- levels(phenos[,phenoCol])

    # Subset
    vcf <- subsetByFilter(x = vcf, filter = filter)

    ggDataList <- bplapply(
        X = pLevels,
        FUN = .densityVepInPhenoLevel,
        vcf = vcf,
        phenoCol = phenoCol,
        vepCol = vepCol,
        param = param,
        filter = VcfFilterRules(), # skip: already filtered above
        unique = unique,
        facet = facet,
        plot = FALSE, # skip: only collect data for each level, here
        popFreq = popFreq,
        BPPARAM = bp(param))

    ggData <- do.call(rbind, ggDataList)

    if (plot){

        ggPlot <- ggplot(
            data = ggData,
            mapping = aes_string(vepCol, colour = phenoCol)) +
            geom_density()

        if (!is.null(facet)){
            ggPlot <- ggPlot + facet_wrap(facets = facet)
        }

        return(ggPlot)

    } else {
        return(ggData)
    }
}

.densityVepInPhenoLevel <- function(
    level, vcf, phenoCol, vepCol, param, filter = VcfFilterRules(),
    unique = FALSE, facet = NULL, plot = FALSE, popFreq = FALSE){

    # Pass relevant namespaces to the parallel environment
    requireNamespace("Biostrings")
    requireNamespace("SummarizedExperiment")
    requireNamespace("IRanges")

    phenos <- colData(vcf)
    stopifnot(level %in% phenos[,phenoCol])

    # Subset

    vcf <- subsetByFilter(x = vcf, filter = filter)

    # Keep the desired consequence for those variants
    ggData <- .vepInPhenoLevel(
        vcf = vcf, phenoCol = phenoCol, level = level, vepCol = vepCol,
        param = param, unique = unique, facet = facet)

    if (nrow(ggData) > 0){
        # add column stating phenotype level
        ggData[,phenoCol] <- level
        # If dealing with population frequency, trim the ALT: prefix
        if (popFreq) {
            ggData[,vepCol] <- gsub(".*:", "", ggData[,vepCol])
        }
        # Coerce to numeric
        ggData[,vepCol] <- as.numeric(ggData[,vepCol])
    }

    if (plot){

        if (nrow(ggData) == 0)
            stop("No data to plot")

        ggPlot <- ggplot(
            data = ggData,
            mapping = aes_string(vepCol, colour = phenoCol)) +
            geom_density()

        if (!is.null(facet)){
            ggPlot <- ggPlot + facet_wrap(facets = facet)
        }

        return(ggPlot)

    } else {
        return(ggData)
    }

}
