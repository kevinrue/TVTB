
# x=ExpandedVCF ----

setMethod(
    "countGenos", c("ExpandedVCF"),
    function(x, genos, pheno = NULL, level = NULL){

        # If a phenotype level is supplied
        if (!is.null(level)){
            # The corresponding phenotype must be supplied
            if (is.null(pheno)){
                stop("pheno required if level supplied")}
            else {
                # Check valid inputs
                stopifnot(.checkPhenoLevel(x, pheno, level))
                # Subset to the samples associated with the phenotype level
                matrixGenos <- geno(x)[["GT"]][,colData(x)[,pheno] == level]
            }
        } else {
            matrixGenos <- geno(x)[["GT"]]
        }

        # Call internal function
        return(.countGenos(matrixGenos, genos))
    }
)

# Checks ----

# Check validity of phenotype level supplied
.checkPhenoLevel <- function(x, pheno, level){

    stopifnot(pheno %in% colnames(colData(x)))

    stopifnot(level %in% colData(x)[,pheno])

    return(TRUE)
}

# Main ----

# x = matrix
# genos = character(n)
.countGenos <- function(x, genos){
    return(as.integer(rowSums(matrix(
        x %in% genos, nrow(x), ncol(x))
    )))
}
