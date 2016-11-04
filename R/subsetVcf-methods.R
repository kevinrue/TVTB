#### VcfFixedFilter ----

setMethod(
    f = "subsetVcf",
    signature = c(x="ExpandedVCF", filter="VcfFixedFilter", param="missing"),
    definition = function(
        x, filter){

        if (!name(filter) %in% colnames(fixed(x)))
            stop(name(filter), " not in fixed(x)")

        if (is.character(value(filter))){
            value(filter) <- sprintf("\"%s\"", value(filter))
        }

        testCmd <- parse(text = sprintf(
            "fixed(x)[,\"%s\"] %s %s",
            name(filter),
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

        if (is.character(value(filter))){
            value(filter) <- sprintf("\"%s\"", value(filter))
        }

        testCmd <- parse(text = sprintf(
            "info(x)[,\"%s\"] %s %s",
            name(filter),
            condition(filter),
            deparse(value(filter))))

        keep <- which(eval(testCmd))

        return(x[keep])
    }
)

#### VcfVepFilter ----

setMethod(
    f = "subsetVcf",
    signature = c(x="ExpandedVCF", filter="VcfVepFilter", param="tSVEParam"),
    definition = function(
        x, filter, param, ..., vep = FALSE){

        return(.subsetVcfVepFilter(vcf = x, filter = filter, param = param))
    }
)

setMethod(
    f = "subsetVcf",
    signature = c(x="ExpandedVCF", filter="VcfVepFilter", param="missing"),
    definition = function(
        x, filter,
        param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
        ..., vep = FALSE){

        # Only required for vep field
        param <- tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA")
        # override defaults (vep)
        param <- .override.tSVEParam(param = param, ...)

        return(.subsetVcfVepFilter(
            vcf = x, filter = filter, param = param, vep = vep))
    }
)

# Subset variants with 1+ predictions passing filter
.subsetVcfVepFilter <- function(vcf, filter, param, vep = FALSE){
    csq <- parseCSQToGRanges(
        x = vcf, VCFRowID = rownames(vcf), info.key = vep(param))

    stopifnot(is.logical(vep))

    # TODO: optionally subset predictions to those passing filters
    # if (vep){
    #     csq <- .subsetVep(csq, filter)
    #     info(vcf)[,vep(param)] <- .CSQfromGRanges(csq)
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

    return(vcf[keep])
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
