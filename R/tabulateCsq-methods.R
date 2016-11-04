#### ByPhenotype() ----

## param = tSVEParam
setMethod(
    f = "tabulateCsqByPhenotype",
    signature = c(vcf="ExpandedVCF", param="tSVEParam"),
    definition = function(
        vcf, phenoCol, csqCol, param, ...,
        unique = FALSE, facet = NULL, plot = FALSE, percentage = FALSE){

        param <- .override.tSVEParam(param, ...)

        .tabulateCsqByPhenotype(
            vcf = vcf, phenoCol = phenoCol, csqCol = csqCol, param = param,
            unique = unique, facet = facet, plot = plot,
            percentage = percentage)
    }
)

## param = missing
setMethod(
    f = "tabulateCsqByPhenotype",
    signature = c(vcf="ExpandedVCF", param="missing"),
    definition = function(
        vcf, phenoCol, csqCol, alts, param = NULL, ...,
        unique = FALSE, facet = NULL, plot = FALSE, percentage = FALSE){

        if (length(alts) < 2)
            stop(
                "length(alts) must be >= 2: ",
                "Heterozygote and Homozygote alternate genotypes")
        # ref will not be used
        param <- tSVEParam(genos = list("", alts[1], alts[2:length(alts)]))
        # Additional tSVEParam overriden
        param <- .override.tSVEParam(param, ...)

        .tabulateCsqByPhenotype(
            vcf = vcf, phenoCol = phenoCol, csqCol = csqCol, param = param,
            unique = unique, facet = facet, plot = plot,
            percentage = percentage)
    }
)

### InPhenoLevel() ----

## param = tSVEParam
setMethod(
    f = "tabulateCsqInPhenoLevel",
    signature = c(vcf="ExpandedVCF", param="tSVEParam"),
    definition = function(
        level, vcf, phenoCol, csqCol, param, ...,
        unique = FALSE, facet = NULL, plot = FALSE, percentage = FALSE){

        param <- .override.tSVEParam(param, ...)

        .tabulateCsqInPhenoLevel(
            level = level, vcf = vcf, phenoCol = phenoCol, csqCol = csqCol,
            param, unique = unique, facet = facet, plot = plot,
            percentage = percentage)
    }
)

## param = tSVEParam
setMethod(
    f = "tabulateCsqInPhenoLevel",
    signature = c(vcf="ExpandedVCF", param="missing"),
    definition = function(
        level, vcf, phenoCol, csqCol, alts, param = NULL, ...,
        unique = FALSE, facet = NULL, plot = FALSE, percentage = FALSE){

        if (length(alts) < 2)
            stop(
                "length(alts) must be >= 2: ",
                "Heterozygote and Homozygote alternate genotypes")
        # ref will not be used
        param <- tSVEParam(genos = list("", alts[1], alts[2:length(alts)]))
        # Additional tSVEParam overriden
        param <- .override.tSVEParam(param, ...)

        .tabulateCsqInPhenoLevel(
            level = level, vcf = vcf, phenoCol = phenoCol, csqCol = csqCol,
            param, unique = unique, facet = facet, plot = plot,
            percentage = percentage)
    }
)

# Private methods ----

.tabulateCsqByPhenotype <- function(
    vcf, phenoCol, csqCol, param,
    unique = FALSE, facet = NULL, plot = FALSE, percentage = FALSE){

    phenos <- colData(vcf)
    pLevels <- levels(phenos[,phenoCol])

    ggDataList <- bplapply(
        X = pLevels,
        FUN = .csqInPhenoLevel,
        vcf = vcf,
        phenoCol = phenoCol,
        csqCol = csqCol,
        param = param,
        unique = unique,
        facet = facet,
        BPPARAM = bp(param)
    )

    # Add column stating pheno level in each data set
    ggDataList <- bpmapply(
        FUN = function(level, df, phenoCol){
            if (nrow(df) > 0){
                df[,phenoCol] <- level
            }
            return(df)
        },
        level = pLevels,
        df = ggDataList,
        MoreArgs = list(phenoCol = phenoCol),
        BPPARAM = bp(param),
        SIMPLIFY = FALSE)

    # Combine into a single data frame
    ggData <- do.call(rbind, ggDataList)

    # Set all factor levels in the combined data set
    if (nrow(ggData) > 0){
        ggData[,phenoCol] <- factor(
            x = ggData[,phenoCol], levels = pLevels)
    }

    if (plot){

        ggPlot <- ggplot(
            data = ggData,
            mapping = aes_string(phenoCol, fill = csqCol))

        if (!is.null(facet)){
            ggPlot <- ggPlot + facet_wrap(facets = facet)
        }
        if (percentage){
            ggPlot <- ggPlot + geom_bar(position = "fill")
        } else {
            ggPlot <- ggPlot + geom_bar()
        }
        return(ggPlot)

    } else {
        longData <- as.data.frame(table(ggData))

        if (is.null(facet))
            dcastRow <- csqCol
        else
            dcastRow <- paste(csqCol, facet, sep = " + ")

        f <- paste(dcastRow, "~", phenoCol)

        wideData <- reshape2::dcast(
            data = longData,
            formula = as.formula(f),
            value.var = "Freq")

        return(wideData)
    }

}

.tabulateCsqInPhenoLevel <- function(
    level, vcf, phenoCol, csqCol, param,
    unique = FALSE, facet = NULL, plot = FALSE, percentage = FALSE){

    # Pass relevant namespaces to the parallel environment
    requireNamespace("Biostrings")
    requireNamespace("SummarizedExperiment")
    requireNamespace("IRanges")

    phenos <- colData(vcf)
    stopifnot(level %in% phenos[,phenoCol])

    ggData <- .csqInPhenoLevel(
        vcf = vcf, phenoCol = phenoCol, level = level, csqCol = csqCol,
        param = param, unique = unique, facet = facet)

    # If data.frame not empty, create a column "Phenotype" stating pLevel
    if (nrow(ggData) > 0){
        ggData[,phenoCol] <- level
    }

    if (plot){

        if (nrow(ggData) == 0)
            stop("No data to plot")

        ggPlot <- ggplot(
            data = ggData,
            mapping = aes_string(phenoCol, fill = csqCol))

        if (!is.null(facet)){
            ggPlot <- ggPlot + facet_wrap(facets = facet)
        }
        if (percentage){
            ggPlot <- ggPlot +
                geom_bar(position = "fill") +
                ylab("Proportion")
        } else {
            ggPlot <- ggPlot + geom_bar()
        }

        return(ggPlot)
    }

    longData <- as.data.frame(table(ggData))

    if (is.null(facet))
        dcastRow <- csqCol
    else
        dcastRow <- paste(csqCol, facet, sep = " + ")

    f <- paste(dcastRow, "~", phenoCol)

    wideData <- dcast(
        data = longData,
        formula = as.formula(f),
        value.var = "Freq")

    return(wideData)
}
