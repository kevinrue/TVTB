
# ByPhenotype() ----

# vcf=ExpandedVCF,paramTVTBparam ----

setMethod(
    "tabulateVepByPhenotype", c("ExpandedVCF", "TVTBparam"),
    function(
        vcf, phenoCol, vepCol, param, ..., filter = VcfFilterRules(),
        unique = FALSE, facet = NULL, plot = FALSE, percentage = FALSE){

        param <- .override.TVTBparam(param, ...)

        .tabulateVepByPhenotype(
            vcf, phenoCol, vepCol, param,
            filter, unique, facet, plot, percentage)
    }
)

# InPhenoLevel() ----

# vcf=ExpandedVCF,paramTVTBparam ----

setMethod(
    "tabulateVepInPhenoLevel", c("ExpandedVCF", "TVTBparam"),
    function(
        level, vcf, phenoCol, vepCol, param, ..., filter = VcfFilterRules(),
        unique = FALSE, facet = NULL, plot = FALSE, percentage = FALSE){

        param <- .override.TVTBparam(param, ...)

        .tabulateVepInPhenoLevel(
            level, vcf, phenoCol, vepCol, param,
            filter, unique, facet, plot, percentage)
    }
)

# Main methods ----

# .tabulateVepByPhenotype ----

# vcf = ExpandedVCF
# phenoCol = character(1)
# vepCol = character(1)
# param = TVTBparam
# filter = <VCF filter rules>
# unique = logical(1)
# facet = character(1) or NULL
# plot = logical(1)
# percentage = logical(1)
.tabulateVepByPhenotype <- function(
    vcf, phenoCol, vepCol, param, filter,
    unique, facet, plot, percentage){

    # Shortcuts
    phenos <- colData(vcf)
    pLevels <- levels(phenos[,phenoCol])

    # Subset
    vcf <- subsetByFilter(vcf, filter)

    # Calculate data in each phenotype level
    ggDataList <- bplapply(
        pLevels,
        .vepInPhenoLevel,
        vcf = vcf,
        phenoCol = phenoCol,
        vepCol = vepCol,
        param = param,
        unique = unique,
        facet = facet,
        BPPARAM = bp(param)
    )

    # Add phenotype level column in each data set
    ggDataList <- bpmapply(
        function(level, df, phenoCol){
            if (nrow(df) > 0){
                df[,phenoCol] <- level
            }
            return(df)
        },
        level = pLevels,
        df = ggDataList,
        MoreArgs = list(phenoCol = phenoCol),
        SIMPLIFY = FALSE,
        BPPARAM = bp(param)
    )

    # Combine data into a single data.frame
    ggData <- do.call(rbind, ggDataList)

    # Set factor levels in the combined data set
    if (nrow(ggData) > 0){
        # Phenotype levels
        ggData[,phenoCol] <- factor(ggData[,phenoCol], pLevels)
    }

    if (plot){
        # Build ggplot (aes, layers, facets)
        ggPlot <- ggplot(ggData, aes_string(phenoCol, fill = vepCol)) +
            scale_x_discrete(drop = FALSE)

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
        # Reshape data in wide format
        longData <- as.data.frame(table(ggData))

        # Prepare the LHS of the formula with/out facet
        dcastRow <- paste(c(vepCol, facet), sep = " + ")

        f <- paste(dcastRow, "~", phenoCol)

        wideData <- dcast(longData, as.formula(f), value.var = "Freq")

        return(wideData)
    }

}

# level = character(1)
# vcf = ExpandedVCF
# phenoCol = character(1)
# vepCol = character(1)
# param = TVTBparam
# filter = <VCF filter rules>
# unique = logical(1)
# facet = character(1) or NULL
# plot = logical(1)
# percentage = logical(1)
.tabulateVepInPhenoLevel <- function(
    level, vcf, phenoCol, vepCol, param, filter,
    unique, facet, plot, percentage){

    # Validate relevant input
    stopifnot(is.character(level))
    stopifnot(length(level) == 1)

    stopifnot(is.character(phenoCol))
    stopifnot(length(phenoCol) == 1)

    stopifnot(is.character(vepCol))
    stopifnot(length(vepCol) == 1)

    stopifnot(is.logical(unique))
    stopifnot(length(unique) == 1)

    stopifnot(is.logical(plot))
    stopifnot(length(plot) == 1)

    stopifnot(is.logical(percentage))
    stopifnot(length(percentage) == 1)

    # Pass relevant namespaces to the parallel environment
    requireNamespace("Biostrings")
    requireNamespace("SummarizedExperiment")
    requireNamespace("IRanges")

    # Shortcut
    phenos <- colData(vcf)
    stopifnot(level %in% phenos[,phenoCol])

    # Subset
    vcf <- subsetByFilter(vcf, filter)

    # Fetch the desired VEP prediction
    ggData <- .vepInPhenoLevel(
        vcf, phenoCol, level, vepCol, param,
        unique = unique, facet = facet)

    if (nrow(ggData) > 0){
        # Add phenotype column
        ggData[,phenoCol] <- level
    }

    if (plot){
        # Build ggplot (aes, layers, facets)
        ggPlot <- ggplot(ggData, aes_string(phenoCol, fill = vepCol)) +
            scale_x_discrete(drop = FALSE)

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
    } else {
        # Reshape data in wide format
        longData <- as.data.frame(table(ggData))

        # Prepare the LHS of the formula with/out facet
        dcastRow <- paste(c(vepCol, facet), sep = " + ")

        f <- paste(dcastRow, "~", phenoCol)

        wideData <- dcast(longData, as.formula(f), value.var = "Freq")

        return(wideData)
    }

}
