.override.tSVEParam <- function(param, ...){
    dots <- list(...)
    n <- names(dots)

    if ("ref" %in% n)
        hRef(param)[[1]] <- dots[["ref"]]

    if ("het" %in% n)
        het(param)[[1]] <- dots[["het"]]

    if ("alt" %in% n)
        hAlt(param)[[1]] <- dots[["alt"]]

    if ("ranges" %in% n)
        ranges(param) <- dots[["ranges"]]

    if ("aaf" %in% n)
        aaf(param) <- dots[["aaf"]]

    if ("maf" %in% n)
        maf(param) <- dots[["maf"]]

    if ("vep" %in% n)
        vep(param) <- dots[["vep"]]

    if ("bp" %in% n)
        bp(param) <- dots[["bp"]]

    validObject(param)
    return(param)
}

# Constructors ----

setMethod(
    f = "initialize",
    signature = c("tSVEParam"),
    definition = function(
        .Object, genos,
        ranges = GRanges(),
        aaf = "AAF", maf = "MAF", vep = "CSQ", bp = SerialParam()){

        # hRef, het, and hAlt genos required
        .checkGenos(x = genos, length = 3)

        # Fill slots with data
        .Object@genos <- genos
        .Object@ranges <- ranges
        .Object@aaf <- aaf
        .Object@maf <- maf
        .Object@vep <- vep
        .Object@bp <- bp

        validObject(.Object)

        return(.Object)
    })

# genos = list
setMethod(
    f = "tSVEParam",
    signature = c(genos="list"),
    definition = function(
        genos,
        ranges = GRanges(),
        aaf = "AAF", maf = "MAF", vep = "CSQ", bp = SerialParam()){

        if (is.null(names(genos))){
            names(genos) <- c("REF", "HET", "ALT")
        } else {
            if (any(names(genos) == "")){
                stop("All elements of genos must be named, or none")
            }
        }

        .tSVEParam(
            genos = genos, ranges = ranges,
            aaf = aaf, maf = maf, vep = vep, bp = bp)
    }
)

# genos = missing
setMethod(
    f = "tSVEParam",
    signature = c(genos="missing"),
    definition = function(
        ref, het, alt,
        ranges = GRanges(),
        aaf = "AAF", maf = "MAF", vep = "CSQ", bp = SerialParam()){

        genos <- list(REF = ref, HET = het, ALT = alt)

        .checkGenos(x = genos, length = 3)

        .tSVEParam(
            genos = genos, ranges = ranges,
            aaf = aaf, maf = maf, vep = vep, bp = bp)
    }
)

# Main method ----

.tSVEParam <- function(
    genos,
    ranges = GRanges(),
    aaf = "AAF", maf = "MAF", vep = "CSQ", bp = SerialParam()){
    new(
        Class = "tSVEParam",
        genos = genos, ranges = ranges,
        aaf = aaf, maf = maf, vep = vep, bp = bp)
}

# Helpers ----

.checkGenos <- function(x, length){
    # For genos(): length 3 required
    # For carrier(): length 2 required
    # For hRef(), het(), hAlt(): length 1 required
    if (length(x) != length)
        stop("length(genos) must equal 3")

    # Must be all named, or none
    if (is.null(names(x)))
        names(x) <- c("REF", "HET", "ALT")
    else {
        if (any(names(x) == ""))
            stop("All elements of genos must be named, or none")
    }

    # hRef, het, hAlt must be character vectors
    if (any(sapply(x, class) != "character"))
        stop("All genos values must be character vectors")

    # cannot overlap
    if (sum(lengths(x)) != length(unique(unlist(x))))
        stop("genos values must not overlap")
}

# Getters and Setters ----

### genos
setMethod(
    f = "genos",
    signature = c("tSVEParam"),
    definition = function(x)
        slot(x, "genos")
)

setReplaceMethod(
    f = "genos", c("tSVEParam", "list"),
    function(x, value){

        .checkGenos(x = value, length = 3)

        slot(x, "genos") <- value
        x
    }
)

### genos
setMethod(
    f = "ranges",
    signature = c("tSVEParam"),
    definition = function(x)
        slot(x, "ranges")
)

setReplaceMethod(
    f = "ranges", c("tSVEParam", "GRanges"),
    function(x, value){
        slot(x, "ranges") <- value
        x
    }
)

### aaf
setMethod(
    f = "aaf",
    signature = c("tSVEParam"),
    definition = function(x)
        slot(x, "aaf")
)

setReplaceMethod(
    f = "aaf", c("tSVEParam", "character"),
    function(x, value){
        if (length(x) != 1)
            stop("length(x) must equal 1")
        slot(x, "aaf") <- value
        x
    }
)

### maf
setMethod(
    f = "maf",
    signature = c("tSVEParam"),
    definition = function(x)
        slot(x, "maf")
)

setReplaceMethod(
    f = "maf", c("tSVEParam", "character"),
    function(x, value){
        if (length(x) != 1)
            stop("length(x) must equal 1")
        slot(x, "maf") <- value
        x
    }
)

### vep
setMethod(
    f = "vep",
    signature = c("tSVEParam"),
    definition = function(x)
        c(slot(x, "vep"))
)

setReplaceMethod(
    f = "vep", c("tSVEParam", "character"),
    function(x, value){
        if (length(value) != 1)
            stop("length(value) must equal 1")
        slot(x, "vep") <- value
        x
    }
)

### hRef
setMethod(
    f = "hRef",
    signature = c("tSVEParam"),
    definition = function(x)
        slot(x, "genos")[1]
)

setReplaceMethod(
    f = "hRef", c("tSVEParam", "list"),
    function(x, value){

        .checkGenos(x = value, length = 1)

        # For unnamed elements, use current value
        if (names(value) == "")
            names(value) <- names(slot(x, "hRef"))

        slot(x, "genos")[1] <- value
        x
    }
)

setReplaceMethod(
    f = "hRef", c("tSVEParam", "character"),
    function(x, value){
        slot(x, "genos")[[1]] <- value
        x
    }
)

### het
setMethod(
    f = "het",
    signature = c("tSVEParam"),
    definition = function(x)
        slot(x, "genos")[2]
)

setReplaceMethod(
    f = "het", c("tSVEParam", "list"),
    function(x, value){

        .checkGenos(x = value, length = 1)

        # For unnamed elements, use current value
        if (names(value) == "")
            names(value) <- names(slot(x, "het"))

        slot(x, "genos")[2] <- value
        x
    }
)

setReplaceMethod(
    f = "het", c("tSVEParam", "character"),
    function(x, value){
        slot(x, "genos")[[2]] <- value
        x
    }
)

### hAlt
setMethod(
    f = "hAlt",
    signature = c("tSVEParam"),
    definition = function(x)
        slot(x, "genos")[3]
)

setReplaceMethod(
    f = "hAlt", c("tSVEParam", "list"),
    function(x, value){

        .checkGenos(x = value, length = 1)

        # For unnamed elements, use current value
        if (names(value) == "")
            names(value) <- names(slot(x, "hAlt"))

        slot(x, "genos")[3] <- value
        x
    }
)

setReplaceMethod(
    f = "hAlt", c("tSVEParam", "character"),
    function(x, value){
        slot(x, "genos")[[3]] <- value
        x
    }
)

### carrier
setMethod(
    f = "carrier",
    signature = c("tSVEParam"),
    definition = function(x)
        c(slot(x, "genos")[2:3])
)

setReplaceMethod(
    f = "carrier", c("tSVEParam", "list"),
    function(x, value){

        .checkGenos(x = value, length = 2)

        slot(x, "genos")[2:3] <- value
        x
    }
)

### bp (BiocParallel)
setMethod(
    f = "bp",
    signature = c("tSVEParam"),
    definition = function(x)
        c(slot(x, "bp"))
)

setReplaceMethod(
    f = "bp", c("tSVEParam", "BiocParallelParam"),
    function(x, value){

        slot(x, "bp") <- value
        x
    }
)
