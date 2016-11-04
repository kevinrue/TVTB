
# tSVEParam ----

.valid.tSVEParam <- function(object){
    g <- slot(object, "genos")

    errors <- c()
    if (length(g) != 3)
        errors <- c(errors, "length(genos) must equal 3")

    if (is.null(names(g)))
        errors <- c(errors, "names(genos) must not be NULL")

    if (sum(lengths(g)) != length(unique(unlist(g))))
        errors <- c(errors, "genos values must not overlap")

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


tSVEParam <- setClass(
    Class = "tSVEParam",

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

    validity = .valid.tSVEParam

)

# VcfBasicFilter ----

.valid.VcfBasicFilter <- function(object){

    errors <- c()

    if (length(slot(object, "name")) != 1)
        errors <- c(errors, "name(x) must equal 3")

    if (length(slot(object, "condition")) != 1)
        errors <- c(errors, "condition(x) must equal 3")

    if (length(errors) > 0)
        return(errors)

    return(TRUE)
}

VcfBasicFilter <- setClass(
    Class = "VcfBasicFilter",

    contains = "VIRTUAL",

    # Define the slots
    slots = c(
        name = "character",
        condition = "character",
        value = "ANY",
        type = "character"
    ),

    validity = .valid.VcfBasicFilter

)

VcfFixedFilter <- setClass(
    Class = "VcfFixedFilter",

    contains = "VcfBasicFilter",

    prototype = list(type = "fixed")

)

VcfInfoFilter <- setClass(
    Class = "VcfInfoFilter",

    contains = "VcfBasicFilter",

    prototype = list(type = "info")

)

VcfVepFilter <- setClass(
    Class = "VcfVepFilter",

    contains = "VcfBasicFilter",

    prototype = list(type = "VEP")

)

# VcfFilterList ----

.valid.VcfFilterList <- function(object){

    validFilterTypes <- c(
        "VcfFixedFilter",
        "VcfInfoFilter",
        "VcfVepFilter")

    errors <- c()

    idxMatch <- match(
        x = sapply(X = filterRules(object), FUN = "class"),
        table = validFilterTypes)
    idxInvalid <- which(is.na(idxMatch))

    if (length(idxInvalid) > 0)
        errors <- c(
            errors,
            paste(
                "invalid filter class for filter(s): ",
                paste(idxInvalid, sep = ", ")
                )
            )

    if (length(errors) > 0)
        return(errors)

    return(TRUE)
}

VcfFilterList <- setClass(
    Class = "VcfFilterList",

    # Define the slots
    slots = c(
        filterRules = "list"
    ),

    validity = .valid.VcfFilterList

)
