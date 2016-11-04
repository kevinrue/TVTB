#### ByPhenotype() ----

## param = tSVEParam
setMethod(
    f = "densityCsqByPhenotype",
    signature = c(vcf="ExpandedVCF", param="tSVEParam"),
    definition = function(
        vcf, phenoCol, csqCol, param, ...,
        unique = FALSE, facet = NULL, plot = FALSE, popFreq = FALSE){

        param <- .override.tSVEParam(param, ...)

        .densityCsqByPhenotype(
            vcf = vcf, phenoCol = phenoCol, csqCol = csqCol, param = param,
            unique = unique, facet = facet, plot = plot, popFreq = popFreq)
    }
)

## param = missing
setMethod(
    f = "densityCsqByPhenotype",
    signature = c(vcf="ExpandedVCF", param="missing"),
    definition = function(
        vcf, phenoCol, csqCol, alts, param = NULL, ...,
        unique = FALSE, facet = NULL, plot = FALSE, popFreq = FALSE){

        if (length(alts) < 2)
            stop(
                "length(alts) must be >= 2: ",
                "Heterozygote and Homozygote alternate genotypes")
        # ref will not be used
        param <- tSVEParam(genos = list("", alts[1], alts[2:length(alts)]))
        # Additional parameters overriden
        param <- .override.tSVEParam(param, ...)

        .densityCsqByPhenotype(
            vcf = vcf, phenoCol = phenoCol, csqCol = csqCol, param = param,
            unique = unique, facet = facet, plot = plot, popFreq = popFreq)
    }
)

#### InPhenoLevel() ----

## param = tSVEParam
setMethod(
    f = "densityCsqInPhenoLevel",
    signature = c(vcf="ExpandedVCF", param="tSVEParam"),
    definition = function(
        level, vcf, phenoCol, csqCol, param, ...,
        unique = FALSE, facet = NULL, plot = FALSE, popFreq = FALSE){

        param <- .override.tSVEParam(param, ...)

        .densityCsqInPhenoLevel(
            level = level, vcf = vcf, phenoCol = phenoCol, csqCol = csqCol,
            param = param,
            unique = unique, facet = facet, plot = plot, popFreq = popFreq)
    }
)

## param = missing
setMethod(
    f = "densityCsqInPhenoLevel",
    signature = c(vcf="ExpandedVCF", param="missing"),
    definition = function(
        level, vcf, phenoCol, csqCol, alts, param = NULL, ...,
        unique = FALSE, facet = NULL, plot = FALSE, popFreq = FALSE){

        if (length(alts) < 2)
            stop(
                "length(alts) must be >= 2: ",
                "Heterozygote and Homozygote alternate genotypes")
        # ref will not be used
        param <- tSVEParam(genos = list("", alts[1], alts[2:length(alts)]))
        # Additional parameters overriden
        param <- .override.tSVEParam(param, ...)

        .densityCsqInPhenoLevel(
            level = level, vcf = vcf, phenoCol = phenoCol, csqCol = csqCol,
            param = param,
            unique = unique, facet = facet, plot = plot, popFreq = popFreq)
    }
)

# Private methods ----

.densityCsqByPhenotype <- function(
    vcf, phenoCol, csqCol, param,
    unique = FALSE, facet = NULL, plot = FALSE, popFreq = FALSE){

    phenos <- colData(vcf)
    pLevels <- levels(phenos[,phenoCol])

    ggDataList <- bplapply(
        X = pLevels,
        FUN = .densityCsqInPhenoLevel,
        vcf = vcf,
        phenoCol = phenoCol,
        csqCol = csqCol,
        param = param,
        unique = unique,
        facet = facet,
        plot = FALSE,
        popFreq = popFreq,
        BPPARAM = bp(param))

    ggData <- do.call(rbind, ggDataList)

    if (plot){

        ggPlot <- ggplot(
            data = ggData,
            mapping = aes_string(csqCol, colour = phenoCol)) +
            geom_density()

        if (!is.null(facet)){
            ggPlot <- ggPlot + facet_wrap(facets = facet)
        }

        return(ggPlot)

    } else {
        return(ggData)
    }
}

.densityCsqInPhenoLevel <- function(
    level, vcf, phenoCol, csqCol, param,
    unique = FALSE, facet = NULL, plot = FALSE, popFreq = FALSE){

    # Pass relevant namespaces to the parallel environment
    requireNamespace("Biostrings")
    requireNamespace("SummarizedExperiment")
    requireNamespace("IRanges")

    phenos <- colData(vcf)
    stopifnot(level %in% phenos[,phenoCol])

    # Keep the desired consequence for those variants
    ggData <- .csqInPhenoLevel(
        vcf = vcf, phenoCol = phenoCol, level = level, csqCol = csqCol,
        param = param, unique = unique, facet = facet)

    if (nrow(ggData) > 0){
        # add column stating phenotype level
        ggData[,phenoCol] <- level
        # If dealing with population frequency, trim the ALT: prefix
        if (popFreq) {
            ggData[,csqCol] <- gsub(".*:", "", ggData[,csqCol])
        }
    }

    # Coerce to numeric
    ggData[,csqCol] <- as.numeric(ggData[,csqCol])

    if (plot){

        if (nrow(ggData) == 0)
            stop("No data to plot")

        ggPlot <- ggplot(
            data = ggData,
            mapping = aes_string(csqCol, colour = phenoCol)) +
            geom_density()

        if (!is.null(facet)){
            ggPlot <- ggPlot + facet_wrap(facets = facet)
        }

        return(ggPlot)

    } else {
        return(ggData)
    }

}
