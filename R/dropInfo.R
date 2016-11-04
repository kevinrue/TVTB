
# main method ----

.dropInfo <- function(vcf, key){

    stopifnot(key %in% rownames(info(header(vcf))))

    idx <- which(colnames(info(vcf)) == key)

    info(vcf) <- info(vcf)[,-idx]
    info(header(vcf)) <- info(header(vcf))[-idx,]

    return(vcf)
}
