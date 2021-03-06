
# vcf=VCF ----

setMethod(
    "plotInfo", c("VCF"),
    function(
        vcf, metric, range, annotation, phenotype,
        type = c("p", "heatmap"), zero.rm = FALSE){
        .plotInfo(vcf, metric, range, annotation, phenotype, type, zero.rm)
    }
)

# Main method ----

.geneRegionTrackFromEnsDb <- function(package, range){
    # Set options (origina value must be restore at the end of the main method)
    options(ucscChromosomeNames = FALSE)
    # TODO: geneRegionTrackFromPkg() S4 method that dispatches EnsDb & TxDb
    # Range-based filter
    edb.filter <- GRangesFilter(range)
    # Fetch all exons by transcript overlapping range
    # TODO: allow showing only gene-level annotations, instead of all exons
    grData <- exonsBy(package, filter = edb.filter)
    # Format data to plot
    unlistData <- unlist(grData)
    unlistData$transcript <- rep(names(grData), lengths(grData))
    unlistData$symbol <- mapIds(
        package, rep(names(grData), lengths(grData)), "GENENAME", "TXID"
    )
    # Build the GeneRegionTrack
    # TODO: allow control of the displayed range (e.g. larger than data)
    grtrack <- GeneRegionTrack(
        # TODO: start = <start of visible region>,
        # TODO: end = <end of visible region>,
        rstarts = start(unlistData),
        rends = end(unlistData),
        chromosome = unique(seqnames(unlistData)),
        transcript = unlistData$transcript,
        gene = unlistData$symbol,
        symbol = unlistData$symbol,
        feature = unlistData$transcript,
        exon = unlist(grData)$exon_id,
        name = "Gene models",
        strand = strand(unlistData),
        showId = TRUE, geneSymbol = TRUE
    )
    return(grtrack)
}

# TODO: allow subset of phenotype levels to be plotted instead of all
.plotInfo <- function(
    vcf, metric, range, annotation, phenotype, type = c("p", "heatmap"),
    zero.rm = FALSE){
    o.ucscChromosomeNames <- getOption("ucscChromosomeNames", FALSE)
    # range must be a GRanges of length 1
    # TODO: allow GRanges of length N, provided all ranges are on 1 sequence
    stopifnot(is(range, "GRanges"))
    stopifnot(length(range) == 1)
    # Identify columns to plot
    metricPattern <- sprintf("^%s_(.*)_%s$", phenotype, metric)
    # TODO: allow phenotype = NA_character_ default to plot a single INFO
    # (metric), instead of one by phenotype level (phenotype_level_metric)
    metricCols <- .findInfoMetricColumns(vcf, metricPattern)
    # Identify levels of phenotype available to plot
    phenoLevels <- gsub(metricPattern, "\\1", metricCols)
    # TODO: Check that all metricCols contain values coercible to numeric
    # Build a GeneRegionTrack from annotation and range
    # TODO: support TxDb
    grTrack <- .geneRegionTrackFromEnsDb(annotation, range)
    # Build a DataTrack from vcf and range
    # TODO: allow choice of plot types
    plotData <- rowRanges(vcf)
    mcols(plotData) <- info(vcf)[,metricCols]
    # Hide zero values if requested
    if (zero.rm){
        # TODO: bplapply
        lData <- lapply(metricCols, function(cName){
            cData <- mcols(plotData)[,cName]
            cData[cData == 0] <- NA
            return(cData)
        })
        for (cIndex in seq_along(metricCols)){
            metricName <- metricCols[cIndex]
            mcols(plotData)[,metricName] <- lData[[cIndex]]
        }
    }
    mafTrack <- DataTrack(
        plotData,
        groups = gsub(metricPattern, "\\1", colnames(mcols(plotData))),
        name = metric,
        type = type
    )
    # Display the plot and return the list value
    p <- plotTracks(
        list(grTrack, mafTrack),
        min(start(ranges(range))),
        max(end(ranges(range)))
    )
    options(ucscChromosomeNames = o.ucscChromosomeNames)
    return(invisible(p))
}
