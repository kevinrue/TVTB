
# Constructors ----

# genos=Genotypes ----

setMethod(
    "TVTBparam", "Genotypes",
    function(
        genos, ranges = GRangesList(),
        aaf = "AAF", maf = "MAF", vep = "CSQ", bp = SerialParam(),
        svp = ScanVcfParam(which = reduce(unlist(ranges)))){

        # Support objects coercible to GRangesList
        ranges <- as(ranges, "GRangesList")

        return(new(
            Class = "TVTBparam",
            genos = genos, ranges = ranges,
            aaf = aaf, maf = maf, vep = vep, bp = bp, svp = svp))
    }
)

# Getters and Setters ----

# genos ----

setMethod("genos", "TVTBparam", function(x) x@genos)

setReplaceMethod(
    "genos", c("TVTBparam", "Genotypes"),
    function(x, value){
        x@genos <- value
        validObject(x)
        return(x)
    }
)

# ranges ----

setMethod("ranges", "TVTBparam", function(x) x@ranges)

setReplaceMethod(
    "ranges", c("TVTBparam", "GRangesList"),
    function(x, value){
        x@ranges <- value
        validObject(x)
        return(x)
    }
)

# aaf ----

setMethod("aaf", "TVTBparam", function(x) x@aaf)

setReplaceMethod(
    "aaf", c("TVTBparam", "character"),
    function(x, value){
        x@aaf <- value
        validObject(x)
        return(x)
    }
)

# maf ----

setMethod("maf", "TVTBparam", function(x) x@maf)

setReplaceMethod(
    "maf", c("TVTBparam", "character"),
    function(x, value){
        x@maf <- value
        validObject(x)
        return(x)
    }
)

# vep ----

setMethod("vep", "TVTBparam", function(x) x@vep)

setReplaceMethod(
    "vep", c("TVTBparam", "character"),
    function(x, value){
        x@vep <- value
        validObject(x)
        return(x)
    }
)

# ref ----

setMethod("ref", "TVTBparam", function(x) ref(x@genos))

setReplaceMethod(
    "ref", c("TVTBparam", "character"),
    function(x, value){
        ref(x@genos) <- value
        validObject(x)
        return(x)
    }
)

# het ----

setMethod("het", "TVTBparam", function(x) het(x@genos))

setReplaceMethod(
    "het", c("TVTBparam", "character"),
    function(x, value){
        het(x@genos) <- value
        validObject(x)
        return(x)
    }
)

# alt ----

setMethod("alt", "TVTBparam", function(x) alt(x@genos))

setReplaceMethod(
    "alt", c("TVTBparam", "character"),
    function(x, value){
        alt(x@genos) <- value
        validObject(x)
        return(x)
    }
)

# carrier ----

setMethod("carrier", "TVTBparam", function(x) carrier(x@genos))

# bp (BiocParallel) ----

setMethod("bp", "TVTBparam", function(x) x@bp)

setReplaceMethod(
    "bp", c("TVTBparam", "BiocParallelParam"),
    function(x, value){
        x@bp <- value
        validObject(x)
        return(x)
    }
)

# suffix ----

setMethod("suffix", "TVTBparam", function(x){
    c(suffix(x@genos), aaf = x@aaf, maf = x@maf)
})

# svp ----

setMethod("svp", "TVTBparam", function(x) x@svp)

setReplaceMethod(
    "svp", c("TVTBparam", "ScanVcfParam"),
    function(x, value){
        x@svp <- value
        validObject(x)
        return(x)
    }
)

# setAs ----

setAs(
    "TVTBparam", "ScanVcfParam",
    function(from){
        return(from@svp)
    }
)
