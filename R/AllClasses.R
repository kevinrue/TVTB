# nocov start
# Genotypes ----

.valid.Genotypes <- function(object){

    errors <- c()

    # Only warn if genotypes overlap
    if (length(genos(object)) != length(unique(genos(object)))){
        warning("genotype values overlap")
        # errors <- c(errors, "genotype values must not overlap")
    }

    if (length(object@ref) == 0){
        warning("No homozygote reference genotypes defined.")
    } else if (any(is.na(object@ref))){
        warning("NA present in homozygote reference genotypes.")
    }

    if (length(object@het) == 0){
        warning("No heterozygote genotypes defined.")
    } else if (any(is.na(object@het))){
        warning("NA present in heterozygote genotypes.")
    }

    if (length(object@alt) == 0){
        warning("No homozygote alternate genotypes defined.")
    } else if (any(is.na(object@alt))){
        warning("NA present in homozygote alternate genotypes.")
    }

    if (any(is.na(match(c("ref", "het", "alt"), names(object@suffix))))){
        errors <- c(
            errors, "names(suffix) must be \"ref\", \"het\", and \"alt\"")
    }

    if (length(errors > 0)){
        return(errors)
    }

    return(TRUE)
}

Genotypes <- setClass(
    "Genotypes",

    # Define the slots
    slots = c(
        ref = "character",
        het = "character",
        alt = "character",
        suffix = "character"
    ),

    validity = .valid.Genotypes
)

# TVTBparam ----

.valid.TVTBparam <- function(object){
    errors <- c()

    if (length(slot(object, "aaf")) != 1){
        errors <- c(errors, "length(aaf) must equal 1")
    }

    if (length(slot(object, "maf")) != 1){
        errors <- c(errors, "length(maf) must equal 1")
    }

    if (length(slot(object, "vep")) != 1){
        errors <- c(errors, "length(vep) must equal 1")
    }

    suffixes <- c(object@genos@suffix, object@aaf, object@maf)
    if (length(suffixes) != length(unique(suffixes))){
        errors <- c(errors, "suffixes must not overlap")
    }

    svpInfo <- vcfInfo(object@svp)
    if ((length(svpInfo) > 0) & all(!object@vep %in% svpInfo)){
        errors <- c(
            errors, "vep %in% vcfWhich(svp) is not TRUE")
    }

    if (length(errors > 0)){
        return(errors)
    }

    return(TRUE)
}


TVTBparam <- setClass(
    "TVTBparam",

    # Define the slots
    slots = c(
        genos = "Genotypes",
        ranges = "GRangesList",
        aaf = "character",
        maf = "character",
        vep = "character",
        bp = "BiocParallelParam",
        svp = "ScanVcfParam"
    ),

    # Set the default values for the slots. (optional)
    prototype = list(
        ranges = GRangesList(),
        aaf = "AAF",
        maf = "MAF",
        vep = "CSQ",
        bp = SerialParam(),
        svp <- ScanVcfParam()
    ),

    validity = .valid.TVTBparam
)

# VcfFilterRules & Co. ----

# So far, no difference from FilterRules
.valid.VcfBasicRules <- function(object){

    errors <- c()

    # TODO: uncomment when relevant validity checks are implemented
    # if (length(errors > 0)){
    #     return(errors)
    # }

    return(TRUE)
}

VcfFixedRules <- setClass(
    "VcfFixedRules",

    contains = "FilterRules",

    validity = .valid.VcfBasicRules
)

VcfInfoRules <- setClass(
    "VcfInfoRules",

    contains = "FilterRules",

    validity = .valid.VcfBasicRules
)

.valid.VcfVepRules <- function(object){
    errors <- c()

    if (!is.character(object@vep)){
        errors <- c(errors, "is.character(vep(x)) is not TRUE")
    }

    if (length(object@vep) != 1){
        errors <- c(errors, "length(vep(x)) must equal 1")
    }

    if (length(errors > 0)){
        return(errors)
    }

    return(TRUE)
}

VcfVepRules <- setClass(
    "VcfVepRules",
    slots = list(vep = "character"),
    prototype = list(vep = "CSQ"),
    contains = "FilterRules",

    validity = .valid.VcfVepRules
)

.valid.VcfFilterRules <- function(object){

    errors <- c()

    validTypes <- c(
        FilterRules = "filter",
        VcfFixedRules = "fixed",
        VcfInfoRules = "info",
        VcfVepRules = "vep")
    typeMatch <- match(object@type, validTypes)

    if (sum(is.na(typeMatch)) > 0){
        errors <- c(errors, "Invalid filter types")
    }

    if (!is.character(object@type)){
        errors <- c(errors, "is.character(type(x)) is not TRUE")
    }

    if (length(object@type) != length(object)){
        errors <- c(errors, "length(type(x)) must equal length(x)")
    }

    if (!is.character(object@vep)){
        errors <- c(errors, "is.character(vep(x)) is not TRUE")
    }

    if (length(object@vep) != 1){
        errors <- c(errors, "length(vep(x)) must equal 1")
    }

    if (length(errors > 0)){
        return(errors)
    }

    return(TRUE)
}

VcfFilterRules <- setClass(
    "VcfFilterRules",
    slots = list(
        type = "character",
        vep = "character"),
    contains = "FilterRules",

    prototype = list(
        type = character(),
        vep = NA_character_
    ),

    validity = .valid.VcfFilterRules
)
# nocov end
