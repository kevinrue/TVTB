
.valid.tSVEParam <- function(object){
    g <- genos(object)

    errors <- c()
    if (length(g) != 3)
        errors <- c(errors, "length(genos) must equal 3")

    if (is.null(names(g)))
        errors <- c(errors, "names(genos) must not be NULL")

    if (any(sapply(g, class) != "character"))
        errors <- c(errors, "All three genos values must be character vectors")

    if (sum(lengths(g)) != length(unique(unlist(g))))
        stop("genos values must not overlap")

    if (length(carrier(object)) != 2)
        errors <- c(errors, "length(carrier) must equal 2")

    if (length(aaf(object)) != 1)
        errors <- c(errors, "length(aaf) must equal 1")

    if (length(maf(object)) != 1)
        errors <- c(errors, "length(maf) must equal 1")

    if (length(vep(object)))
        errors <- c(errors, "length(vep) must equal 1")

    NULL
}


tSVEParam <- setClass(
    Class = "tSVEParam",

    # Define the slots
    slots = c(
        genos = "list",
        aaf = "character",
        maf = "character",
        vep = "character",
        bp = "BiocParallelParam"
        # TODO: GQ, MAF, ...
    ),

    # Set the default values for the slots. (optional)
    prototype = list(
        aaf = "AAF",
        maf = "MAF",
        vep = "CSQ",
        bp = SerialParam()
        # TODO: GQ, MAF, ...
    ),

    validity=.valid.tSVEParam

)

