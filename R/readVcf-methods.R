
.checkParam <- function(param, file, ...){
    vh <- scanVcfHeader(file)

    phenoSamples <- rownames(colData)

    if(!vep(param) %in% rownames(info(vh))){
        stop("\"", vep(param), "\" not found among VCF INFO keys.")
    }

    SVP <- as(param, "ScanVcfParam")

    return(SVP)
}

# Signatures ----

# NOTE: Signature file=ANY,param=TVTBparam causes following error

# Note: method with signature 'character#ANY' chosen for function 'readVcf',
# target signature 'character#TVTBparam'.
# "ANY#TVTBparam" would also be valid
# Error in (function (classes, fdef, mtable)  :
#     unable to find an inherited method for function 'seqinfo' for
#     signature '"TVTBparam"'

setMethod(
    "readVcf", c("character", "TVTBparam"),
    function(file, genome, param, ..., colData = DataFrame()){
        .readVcf(file, genome, param, ..., colData = colData)
    }
)

setMethod(
    "readVcf", c("TabixFile", "TVTBparam"),
    function(file, genome, param, ..., colData = DataFrame()){
        .readVcf(file, genome, param, colData = colData, ...)
    }
)

# Main method ----

.readVcf <- function(file, genome, param, ..., colData = DataFrame()){

    # Check TVTBparam and coerce to ScanVcfParam
    SVP <- .checkParam(param, file, ...)

    # Import from VCF only samples given in phenotypes
    if (nrow(colData) > 0){
        if (length(vcfSamples(SVP)) > 0){
            if (length(vcfSamples(SVP)) != nrow(colData) ||
                any(sort(vcfSamples(SVP)) != sort(rownames(colData)))){
                stop("Different samples in svp(param) and colData.")
            }
        }
        vcfSamples(SVP) <- rownames(colData)
    }

    # Use VariantAnnotation readVcf
    vcf <- readVcf(file, genome, SVP, ...)

    # Attach (default) phenotype information
    if (nrow(colData) > 0){
        colData(vcf) <- colData
    }

    metadata(vcf)[["TVTBparam"]] <- param

    return(vcf)
}
