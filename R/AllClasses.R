
# TVTBparam ----

.valid.TVTBparam <- function(object){
    g <- slot(object, "genos")

    errors <- c()
    if (length(g) != 3)
        errors <- c(errors, "length(genos) must equal 3")

    if (is.null(names(g)))
        errors <- c(errors, "names(genos) must not be NULL")

    if (any(lengths(g) < 1))
        errors <- c(errors, "lengths(genos) must all be >= 1")

    if (sum(lengths(g)) != length(unique(unlist(g))))
        errors <- c(errors, "genos values must not overlap")

    suffixes <- c(names(g), object@aaf, object@maf)
    if (length(suffixes) != length(unique(suffixes)))
        errors <- c(errors, "suffixes must not overlap")

    if (length(slot(object, "aaf")) != 1)
        errors <- c(errors, "length(aaf) must equal 1")

    if (length(slot(object, "maf")) != 1)
        errors <- c(errors, "length(maf) must equal 1")

    if (length(slot(object, "vep")) != 1)
        errors <- c(errors, "length(vep) must equal 1")

    if (length(errors > 0))
        return(errors)

    return(TRUE)
}


TVTBparam <- setClass(
    Class = "TVTBparam",

    # Define the slots
    slots = c(
        genos = "list",
        ranges = "GRanges",
        aaf = "character",
        maf = "character",
        vep = "character",
        bp = "BiocParallelParam"
        # TODO: GQ, MAF, ...
    ),

    # Set the default values for the slots. (optional)
    prototype = list(
        ranges = GRanges(),
        aaf = "AAF",
        maf = "MAF",
        vep = "CSQ",
        bp = SerialParam()
        # TODO: GQ, MAF, ...
    ),

    validity = .valid.TVTBparam

)

# VcfFilterRules & Co. ----

# TODO: any difference from FilterRules?
.valid.VcfBasicRules <- function(object){
    errors <- c()

    rulesName <- names(object@listData)
    if (length(rulesName) != length(unique(rulesName)))
        errors <- c(errors, "names must be unique")

    if (length(errors > 0))
        return(errors)

    return(TRUE)
}

VcfFixedRules <- setClass(
    Class = "VcfFixedRules",

    contains = "FilterRules",

    validity = .valid.VcfBasicRules
)


VcfInfoRules <- setClass(
    Class = "VcfInfoRules",

    contains = "FilterRules",

    validity = .valid.VcfBasicRules
)

.valid.VcfVepRules <- function(object){
    errors <- c()

    if (length(vep(object)) != 1)
        errors <- c(errors, "length(vep(x)) must equal 1")

    if (class(vep(object)) != "character")
        errors <- c(errors, "vep(x) must be \"character\"")

    rulesName <- names(object@listData)
    if (length(rulesName) != length(unique(rulesName)))
        errors <- c(errors, "names must be unique")

    if (length(errors > 0))
        return(errors)

    return(TRUE)
}

VcfVepRules <- setClass(
    Class = "VcfVepRules",
    slots = list(vep = "character"),
    contains = "FilterRules",

    validity = .valid.VcfVepRules
)

.valid.VcfFilterRules <- function(object){

    errors <- c()

    validTypes <- c("fixed", "info", "vep")
    typeMatch <- match(type(object), validTypes)

    if (sum(is.na(typeMatch)) > 0)
        errors <- c(errors, "one or more invalid filter types")

    if (length(type(object)) != length(object))
        errors <- c(errors, "length(type(x)) must equal length(x)")

    rulesName <- names(object@listData)
    if (length(rulesName) != length(unique(rulesName)))
        errors <- c(errors, "names must be unique")

    if (length(vep(object)) != 1)
        errors <- c(errors, "length(vep(x)) must equal 1")

    if (length(errors > 0))
        return(errors)

    return(TRUE)

}

VcfFilterRules <- setClass(
    Class = "VcfFilterRules",
    slots = list(
        type = "character",
        vep = "character"),
    contains = "FilterRules",

    validity = .valid.VcfFilterRules
)
