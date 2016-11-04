# genotypes = ExpandedVCF
setMethod(
    f = "countGenos",
    signature = c(x="ExpandedVCF"),
    definition = function(x, genos, pheno = NULL, level = NULL){
        # If a phenotype level is supplied
        if (!is.null(level)){
            # The corresponding phenotype must be supplied
            if (is.null(pheno)){
                stop("pheno required if level supplied")}
            else {
                # Check valid inputs
                .checkPhenoLevel(
                    x = x,
                    pheno = pheno,
                    level = level
                )
                # Subset to the samples associated with the phenotype level
                matrixGenos <- geno(x)[["GT"]][,colData(x)[,pheno] == level]
            }
        } else {
            matrixGenos <- geno(x)[["GT"]]
        }
        # Call internal function
        .countGenos(x = matrixGenos, genos = genos)
    }
)

# Main ----

.checkPhenoLevel <- function(x, pheno, level){

    if (! pheno %in% colnames(colData(x))){
        stop("pheno not found in colnames(colData(x))")
    }

    if (! level %in% colData(x)[,pheno]){
        stop("level not found in colData(x)[,pheno]")
    }

    return(TRUE)
}

# x: matrix
# genos: character vector
.countGenos <- function(x, genos){
    as.integer(rowSums(matrix(
        x %in% genos, nrow = nrow(x), ncol = ncol(x))))
}
