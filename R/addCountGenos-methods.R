
# samples as default (all): probably the most common usage
setMethod(
    f = "addCountGenos",
    signature = c(vcf="ExpandedVCF", samples="missing"),
    definition = function(
        vcf, genos, key, description, samples = 1:ncol(vcf), force = FALSE){
        .addCountGenos(
            vcf, genos, key, description, samples = 1:ncol(vcf),
            force = force)}
)

# samples as column index
setMethod(
    f = "addCountGenos",
    signature = c(vcf="ExpandedVCF", samples="numeric"),
    definition = function(
        vcf, genos, key, description, samples = 1:ncol(vcf), force = FALSE){
        .addCountGenos(
            vcf, genos, key, description, samples = 1:ncol(vcf),
            force = force)}
)

# samples as column name
setMethod(
    f = "addCountGenos",
    signature = c(vcf="ExpandedVCF", samples="character"),
    definition = function(
        vcf, genos, key, description, samples = 1:ncol(vcf), force = FALSE){
        .addCountGenos(
            vcf, genos, key, description, samples = 1:ncol(vcf),
            force = force)}
)

.addCountGenos <- function(
    vcf, genos, key, description, samples = 1:ncol(vcf), force = FALSE){

    newInfo <- DataFrame(
        Number = 1,
        Type = "Integer",
        Description = description,
        row.names = key)

    keyIndex <- which(rownames(info(header(vcf))) == key)

    # Replace of append field in header
    if (length(keyIndex) > 0){
        if (force){
            message(key, " will be overwritten.")
        } else {
            stop(key, " field already present.")
        }
    } else {
        keyIndex <- nrow(info(header(vcf))) + 1
    }
    info(header(vcf))[keyIndex,] <- newInfo

    # Replace or add values in INFO
    info(vcf)[,key] <- .countGenos(
        x = geno(vcf)[["GT"]][,samples], genos = genos)

    return(vcf)
}
