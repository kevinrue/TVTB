# TODO: rename level -> phenoLevel?

# ByPhenotype() ----

# vcf=ExpandedVCF,paramTVTBparam ----

setMethod(
    "densityVepByPhenotype", c("ExpandedVCF", "TVTBparam"),
    function(
        vcf, phenoCol, vepCol, param, ..., filter = VcfFilterRules(),
        unique = FALSE, facet = NULL, plot = FALSE, pattern = NULL,
        layer = "density+dotplot"){

        # ellipsis may be used to override some parameters
        param <- .override.TVTBparam(param, ...)

        return(.densityVepByPhenotype(
            vcf, phenoCol, vepCol, param,
            filter, unique, facet, plot, pattern, layer))
    }
)

# InPhenoLevel() ----

# vcf=ExpandedVCF,paramTVTBparam ----

setMethod(
    "densityVepInPhenoLevel", c("ExpandedVCF", "TVTBparam"),
    function(
        level, vcf, phenoCol, vepCol, param, ..., filter = VcfFilterRules(),
        unique = FALSE, facet = NULL, plot = FALSE, pattern = NULL,
        layer = "density+dotplot"){

        # ellipsis may be used to override some parameters
        param <- .override.TVTBparam(param, ...)

        return(.densityVepInPhenoLevel(
            level, vcf, phenoCol, vepCol, param,
            filter, unique, facet, plot, pattern, layer))
    }
)

# ggplot output ----

# Build the ggplot output
# ggData = data.frame
# vepCol = character(1)
# phenoCol = character(1)
# layer = character(1)
# facet = character(1) or NULL
.densityVepPlot <- function(ggData, vepCol, phenoCol, layer, facet){

    # Build basic ggplot object
    ggPlot <- ggplot(ggData, aes_string(vepCol, colour = phenoCol)) +
        scale_colour_discrete(phenoCol)

    # Add requested layers (possibility of None!)
    ggLayers <- strsplit(gsub(" ", "", layer), "+", fixed = TRUE)[[1]]
    if (any(grepl("density", ggLayers))){
        ggPlot <- ggPlot + geom_density()
    }
    if (any(grepl("dotplot", ggLayers))){
        ggPlot <- ggPlot +
            geom_dotplot(
                aes_string(fill = phenoCol),
                binpositions = "all",
                stackgroups = TRUE) +
            scale_fill_discrete(phenoCol)
    }

    if (!is.null(facet)){
        ggPlot <- ggPlot + facet_wrap(facets = facet)
    }

    if (!is.null(facet)){
        ggPlot <- ggPlot + facet_wrap(facets = facet)
    }

    # return ggplot object
    return(ggPlot)
}

# Main methods ----

# .densityVepByPhenotype ----

# vcf = ExpandedVCF
# phenoCol = character(1)
# vepCol = character(1)
# param = TVTBparam
# filter = <VCF filter rules>
# unique = logical(1)
# facet = character(1) or NULL
# plot = logical(1)
# pattern = character(1) including at least one pair of "()" or "
# layer = character(1)
.densityVepByPhenotype <- function(
    vcf, phenoCol, vepCol, param,
    filter, unique, facet, plot, pattern, layer){

    # Shortcuts
    phenos <- colData(vcf)
    pLevels <- levels(phenos[,phenoCol])

    # Subset
    vcf <- subsetByFilter(vcf, filter)

    # Calculate data in each phenotype level
    ggDataList <- bplapply(
        pLevels,
        .densityVepInPhenoLevel,
        vcf = vcf,
        phenoCol = phenoCol,
        vepCol = vepCol,
        param = param,
        filter = VcfFilterRules(), # skip: already filtered above
        unique = unique,
        facet = facet,
        plot = FALSE, # skip: only collect data for each level, here
        pattern = pattern,
        BPPARAM = bp(param))

    # Combine data into a single data.frame
    ggData <- do.call(rbind, ggDataList)

    if (plot){
        # Build ggplot (aes, layers, facets)
        ggPlot <- .densityVepPlot(ggData, vepCol, phenoCol, layer, facet)

        return(ggPlot)

    } else {
        # Retun data in long format
        return(ggData)
    }
}

# .densityVepInPhenoLevel ----

# level = character(1)
# vcf = ExpandedVCF
# phenoCol = character(1)
# vepCol = character(1)
# param = TVTBparam
# filter = <VCF filter rules>
# unique = logical(1)
# facet = character(1) or NULL
# plot = logical(1)
# pattern = character(1) including at least one pair of "()" or "
# layer = character(1)
.densityVepInPhenoLevel <- function(
    level, vcf, phenoCol, vepCol, param, filter = VcfFilterRules(),
    unique = FALSE, facet = NULL, plot = FALSE, pattern = NULL,
    layer = "density+dotplot"){

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

    stopifnot(is.character(layer))
    stopifnot(length(layer) == 1)

    # Pass relevant namespaces to the parallel environment
    requireNamespace("Biostrings")
    requireNamespace("SummarizedExperiment")
    requireNamespace("IRanges")

    # Shortcut
    phenos <- colData(vcf)
    stopifnot(level %in% phenos[,phenoCol])

    # Subset
    vcf <- subsetByFilter(x = vcf, filter = filter)

    # Fetch the desired VEP prediction
    ggData <- .vepInPhenoLevel(
        vcf, phenoCol, level, vepCol, param, unique, facet)

    if (nrow(ggData) > 0){
        # Add phenotype column
        ggData[,phenoCol] <- level
        # Extract requested pattern from VEP prediction
        if (!is.null(pattern)) {
            ggData[,vepCol] <- gsub(pattern, "\\1", ggData[,vepCol])
        }
        # Coerce to numeric
        ggData[,vepCol] <- as.numeric(ggData[,vepCol])
    }

    if (plot){
        # Build ggplot object
        ggPlot <- .densityVepPlot(ggData, vepCol, phenoCol, layer, facet)

        return(ggPlot)

    } else {
        # Return data in long format
        return(ggData)
    }

}
