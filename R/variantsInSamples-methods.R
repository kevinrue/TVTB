
## samples = column index ----

setMethod(
    f = "variantsInSamples",
    signature = c("ExpandedVCF", "TVTBparam"),

    definition = function(
        vcf, param, samples = 1:ncol(vcf), ..., unique = FALSE){

        param <- .override.TVTBparam(param, ...)

        .variantsInSamples(
            vcf = vcf, param = param, samples = samples, unique = unique)}
)

# Main method ----

.variantsInSamples <- function(
    vcf, param, samples = 1:ncol(vcf), unique = FALSE){

    stopifnot(class(samples) %in% c("character", "numeric", "integer"))
    stopifnot(is.logical(unique))
    stopifnot(length(samples) > 0)
    
    # If samples is character (colnames)
    if (identical(class(samples), "character")){
        # convert to index
        samples <- match(samples, colnames(vcf))
        # If some sample names were not found, throw an error
        stopifnot(!any(is.na(samples)))
    }
    
    # Extract once genotypes from ExpandedVCF
    genos <- geno(vcf)[["GT"]]

    # Genotype of samples (selected by index/names)
    genosSamples <- genos[,samples]

    # Count how many samples carry the variant
    countGenoWithPheno <- apply(
        genosSamples, 1, function(x){
            sum(x %in% unlist(carrier(param)))})

    # If all other sample have to be non-carrier
    if (unique){
        # Genotype of other samples
        genosOtherSamples <- genos[,-samples]
        # Count how many samples carry the variant
        countGenoWithoutPheno <- apply(
            genosOtherSamples, 1, function(x){
                sum(x %in% unlist(carrier(param))) > 0})
        # Keep only the variants unique to samples (>0 in samples, 0 in others)
        variantsIdx <- which(
            (countGenoWithPheno / countGenoWithoutPheno) == Inf)
    } else {
        # Otherwise, keep all the variants seen in samples
        variantsIdx <- which(countGenoWithPheno > 0)
    }

    return(variantsIdx)
}
