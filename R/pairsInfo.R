
# vcf=VCF ----

setMethod(
    "pairsInfo", c("VCF"),
    function(vcf, metric, phenotype, ..., title = metric){
        .pairsInfo(vcf, metric, phenotype, ..., title = title)
    }
)

# Main method ----

# TODO: allow subset of phenotype levels to be plotted instead of all
.pairsInfo <- function(vcf, metric, phenotype, ..., title = metric){
    # Identify columns to plot
    metricPattern <- sprintf("^%s_(.*)_%s$", phenotype, metric)
    metricCols <- .findInfoMetricColumns(vcf, metricPattern)
    # Identify levels of phenotype available to plot
    phenoLevels <- gsub(metricPattern, "\\1", metricCols)
    # TODO: Check that all metricCols contain values coercible to numeric
    # Build the data frame to plot
    plotData <- as.data.frame(info(vcf)[,metricCols])
    # Simplify the facet names to the phenotype level
    colnames(plotData) <- gsub(metricPattern, "\\1", colnames(plotData))
    # Return the ggpairs object (or plot)
    ggpairs(plotData, title = title, ...)
}
