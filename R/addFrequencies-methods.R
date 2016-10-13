
# vcf=ExpandedVCF,phenos=list ----

setMethod(
    "addFrequencies", c("ExpandedVCF", "list"),
    function(vcf, phenos, force = FALSE){

        return(.addFrequencies(vcf, phenos, force))
    }
)

# vcf=ExpandedVCF,phenos=character ----

setMethod(
    "addFrequencies", c("ExpandedVCF", "character"),
    function(vcf, phenos, force = FALSE){

        # List all levels of phenotypes supplied
        phenos <- sapply(
            phenos, function(x){unique(colData(vcf)[,x])},
            simplify = FALSE)

        return(.addFrequencies(vcf, phenos, force))
    }
)

# vcf=ExpandedVCF,phenos=missing ----

setMethod(
    "addFrequencies", c("ExpandedVCF", "missing"),
    function(vcf, force = FALSE){

        # Delegate to .addOverallFrequencies if no phenotype is supplied
        return(.addOverallFrequencies(vcf, force))
    }
)

# Main method ----

# vcf = ExpandedVCF
# phenos = list()
# force = logical(1)
.addFrequencies <- function(vcf, phenos, force){

    # Validate relevant inputs
    stopifnot(is.list(phenos))

    # TODO: consider parallel processing followed by merging
    # Iteratively add calculated values for each phenotype level
    for (pheno in names(phenos)){
        for (level in phenos[[pheno]]){
            vcf <- .addPhenoLevelFrequencies(vcf, pheno, level, force)
        }
    }

    return(vcf)
}
