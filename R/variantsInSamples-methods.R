
# vcf=ExpandedVCF ----

setMethod(
    "variantsInSamples", c("ExpandedVCF"),
    function(vcf, samples = 1:ncol(vcf), unique = FALSE){
        return(.variantsInSamples(vcf, samples, unique))
    }
)

# Main method ----

# vcf = ExpandedVCF
# samples = <indices or names to identify samples>
# unique = logical(1)
.variantsInSamples <- function(vcf, samples, unique){

    # Check input
    stopifnot(class(samples) %in% c("character", "numeric", "integer"))
    stopifnot(is.logical(unique))
    stopifnot(length(samples) > 0)
    stopifnot(!any(duplicated(samples)))

    stopifnot("TVTBparam" %in% names(metadata(vcf)))
    param <- metadata(vcf)[["TVTBparam"]]
    validObject(param@genos)

    # Convert sample names to indices: required to select the others later
    if (is.character(samples)){
        # Check that all samples are present
        stopifnot(all(samples %in% colnames(vcf)))
        # Convert to index
        samples <- match(samples, colnames(vcf))
    }

    # Extract once genotypes from ExpandedVCF (shortcut)
    genos <- geno(vcf)[["GT"]]

    # Genotype of samples (selected by index/names)
    genosSamples <- genos[,samples]

    # Count how many samples carry the variant
    genoWithPheno <- apply(
        genosSamples, 1, function(x){
            sum(x %in% unlist(carrier(param))) > 0})

    # If all other sample have to be non-carrier (i.e. unique to samples)
    if (unique){
        # Genotype of other samples
        genosOtherSamples <- genos[,-samples]
        # Count how many samples carry the variant
        genoWithoutPheno <- apply(
            genosOtherSamples, 1, function(x){
                sum(x %in% unlist(carrier(param))) > 0})
        print(genoWithPheno)
        print(genoWithoutPheno)
        print(genoWithPheno & !genoWithoutPheno)
        # Keep only the variants unique to samples
        return(which(genoWithPheno & !genoWithoutPheno)
        )
    }
    # Otherwise, keep all the variants seen in samples
    return(which(genoWithPheno))
}
