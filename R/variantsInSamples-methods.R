### param = tSVEParam ----
## samples = column index ----

setMethod(
    f = "variantsInSamples",
    signature = c("ExpandedVCF", "numeric", "tSVEParam"),

    definition = function(vcf, samples, param = NULL, ..., unique = FALSE){

        param <- .override.tSVEParam(param, ...)

        .variantsInSamples(
            vcf = vcf, samples = samples, param = param,
            unique = unique)}
)

# samples as column names
setMethod(
    f = "variantsInSamples",
    signature = c("ExpandedVCF", "character", "tSVEParam"),

    definition = function(vcf, samples, param = NULL, ..., unique = FALSE){

        param <- .override.tSVEParam(param, ...)

        .variantsInSamples(
            vcf = vcf, samples = samples, param = param,
            unique = unique)}
)

### param = missing ----
## samples = column index ----

setMethod(
    f = "variantsInSamples",
    signature = c("ExpandedVCF", "numeric", "missing"),

    definition = function(
        vcf, samples, alts, unique = FALSE){

        if (length(alts) < 2)
            stop(
                "length(alts) must be >= 2: ",
                "Heterozygote and Homozygote alternate genotypes")
        # ref will not be used
        param <- tSVEParam(genos = list("", alts[1], alts[2:length(alts)]))

        .variantsInSamples(
            vcf = vcf, samples = samples, param = param,
            unique = unique)}
)

# samples as column names
setMethod(
    f = "variantsInSamples",
    signature = c("ExpandedVCF", "character", "missing"),

    definition = function(
        vcf, samples, alts, unique = FALSE){

        if (length(alts) < 2)
            stop(
                "length(alts) must be >= 2: ",
                "Heterozygote and Homozygote alternate genotypes")
        # ref will not be used
        param <- tSVEParam(genos = list("", alts[1], alts[2:length(alts)]))

        .variantsInSamples(
            vcf = vcf, samples = samples, param = param,
            unique = unique)}
)

.variantsInSamples <- function(
    vcf, samples, param, unique = FALSE){

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
