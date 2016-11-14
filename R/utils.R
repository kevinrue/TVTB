
.findInfoMetricColumns <- function(vcf, pattern){
    # Check that at least 1 <phenotype>_*_<metric> exist in info(vcf)
    infoCols <- colnames(info(vcf))
    metricCols <- grep(pattern, infoCols, value = TRUE)
    # Check that at least one INFO column is detected
    if (!length(metricCols)){
        stop("No INFO column found matching pattern ", pattern)
    }
    return(metricCols)
}
