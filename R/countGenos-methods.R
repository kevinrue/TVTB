# genotypes = ExpandedVCF
setMethod(
    f = "countGenos",
    signature = c(x="ExpandedVCF", genos="character"),
    definition = function(x, genos, pheno = NULL, level = NULL){
        if (!is.null(level))
            if (is.null(pheno)){
                stop("level argument is invalid if pheno=NULL")}
            else {
                .checkPhenoLevel(x = x, pheno = pheno, level = level)
                # Subset to the samples associated with the phenotype level
                matrixGenos <- geno(x)[["GT"]][,colData(x)[,pheno] == level]
            }

        .countGenos(
            x = matrixGenos, genos = genos, pheno = pheno, level = level)
    }
)

# genotypes = matrix
setMethod(
    f = "countGenos",
    signature = c(x="matrix", genos="character"),
    definition = function(x, genos){
        .countGenos(x, genos, pheno = NULL, level = NULL)
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


.countGenos <- function(x, genos, pheno = NULL, level = NULL){
    as.integer(rowSums(matrix(
        x %in% genos, nrow = nrow(x), ncol = ncol(x))))
}
