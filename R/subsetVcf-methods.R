#### VcfFixedFilter ----

setMethod(
    f = "subsetVcf",
    signature = c(x="ExpandedVCF", filter="VcfFixedFilter", param="missing"),
    definition = function(
        x, filter){

        if (!name(filter) %in% colnames(fixed(x)))
            stop(name(filter), " not in fixed(x)")

        testCmd <- parse(text = sprintf(
            "fixed(x)[,%s] %s %s",
            deparse(name(filter)),
            condition(filter),
            deparse(value(filter))))

        keep <- which(eval(testCmd))

        return(x[keep])
    }
)

#### VcfInfoFilter ----

setMethod(
    f = "subsetVcf",
    signature = c(x="ExpandedVCF", filter="VcfInfoFilter", param="missing"),
    definition = function(
        x, filter){

        if (!name(filter) %in% colnames(info(x)))
            stop(name(filter), " not in info(x)")

        testCmd <- parse(text = sprintf(
            "info(x)[,%s] %s %s",
            deparse(name(filter)),
            condition(filter),
            deparse(value(filter))))

        keep <- which(eval(testCmd))

        return(x[keep])
    }
)

#### VcfVepFilter ----

setMethod(
    f = "subsetVcf",
    signature = c(x="ExpandedVCF", filter="VcfVepFilter", param="TVTBparam"),
    definition = function(
        x, filter, param, ..., vep = FALSE){

        return(.subsetVcfVepFilter(
            x = x, filter = filter, param = param, vep = vep))
    }
)

setMethod(
    f = "subsetVcf",
    signature = c(x="ExpandedVCF", filter="VcfVepFilter", param="missing"),
    definition = function(
        x, filter,
        param = TVTBparam(ref = "refNA", het = "hetNA", alt = "altNA"),
        ..., vep = FALSE){

        # Only required for vep field
        param <- TVTBparam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.TVTBparam(param = param, ...)

        return(.subsetVcfVepFilter(
            x = x, filter = filter, param = param, vep = vep))
    }
)

# Subset variants with 1+ predictions passing filter
.subsetVcfVepFilter <- function(x, filter, param, vep = FALSE){
    csq <- parseCSQToGRanges(
        x = x, VCFRowID = rownames(x), info.key = vep(param))

    # stopifnot(is.logical(vep))
    # TODO: optionally subset predictions to those passing filters
    # Requires parsingCsqToGRanges, but then... GRangesToCsq!
    # if (vep){
    #     csq <- .subsetVep(csq, filter)
    #     info(x)[,vep(param)] <- .CSQfromGRanges(csq)
    # }

    if (!name(filter) %in% colnames(mcols(csq)))
        stop(name(filter), " not in parseCSQToGRanges(x)")

    testCmd <- parse(text = sprintf(
        paste(
            "subset(x = mcols(csq),",
            "subset = %s %s %s,",
            "select = \"VCFRowID\",",
            "drop = TRUE)"
        ),
        name(filter),
        condition(filter),
        deparse(value(filter))))

    keep <- unique(eval(testCmd))

    return(x[keep])
}

#### VcfFilterList ----

setMethod(
    f = "subsetVcf",
    signature = c(x="ExpandedVCF", filter="VcfFilterList", param="TVTBparam"),
    definition = function(
        x, filter, param, ..., vep = FALSE){

        return(.subsetVcfFilterList(
            x = x, filter = filter, param = param, vep = vep))
    }
)

setMethod(
    f = "subsetVcf",
    signature = c(x="ExpandedVCF", filter="VcfFilterList", param="missing"),
    definition = function(
        x, filter,
        param = TVTBparam(ref = "refNA", het = "hetNA", alt = "altNA"),
        ..., vep = FALSE){

        # Only required for vep field
        param <- TVTBparam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.TVTBparam(param = param, ...)

        return(.subsetVcfFilterList(
            x = x, filter = filter, param = param, vep = vep))
    }
)

.subsetVcfFilterList <- function(x, filter, param, vep){

    for (basicFilter in filterRules(filter)){
        x <- switch (class(basicFilter),
            VcfFixedFilter = subsetVcf(
                x = x, filter = basicFilter),
            VcfInfoFilter = subsetVcf(
                x = x, filter = basicFilter),
            VcfVepFilter = subsetVcf(
                x = x, filter = basicFilter, param = param, vep = vep))
    }

    return(x)
}

# Subset prediction passing the filter
# .subsetVep <- function(csq, filter){
#     testCmd <- parse(text = sprintf(
#         paste(
#             "subset(x = mcols(csq),",
#             "subset = %s %s %s,",
#             "select = !grepl(pattern=\"VCFRowID\", x=colnames(mcols(csq))),",
#             "drop = TRUE)"
#         ),
#         name(filter),
#         condition(filter),
#         value(filter)))
#
#     keep <- eval(testCmd)
#
#     return(csq[keep])
# }
