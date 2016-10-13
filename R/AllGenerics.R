
# Genotypes class ----

setGeneric(
    "Genotypes", signature = c("ref", "het", "alt"),
    function(ref, het, alt, suffix = c(ref="REF", het="HET", alt="ALT"))
        standardGeneric("Genotypes")
)

# dropInfo ----

setGeneric(
    "dropInfo", signature = "vcf",
    function(vcf, key = NULL, slot = "both")
        standardGeneric("dropInfo")
)

# addOverallFrequencies ----

setGeneric(
    "addOverallFrequencies", signature = c("vcf"),
    function(vcf, force = FALSE)
        standardGeneric("addOverallFrequencies")
)

# addPhenoLevelFrequencies ----

setGeneric(
    "addPhenoLevelFrequencies", signature = c("vcf"),
    function(vcf, pheno, level, force = FALSE)
        standardGeneric("addPhenoLevelFrequencies")
)

# addPhenotypeFrequencies ----

setGeneric(
    "addFrequencies", signature = c("vcf", "phenos"),
    function(vcf, phenos, force = FALSE)
        standardGeneric("addFrequencies")
)

# addCountGenos ----

setGeneric(
    "addCountGenos", signature = "vcf",
    function(
        vcf, genos, key, description, samples = 1:ncol(vcf), force = FALSE)
        standardGeneric("addCountGenos")
)

# countGenos ----

setGeneric(
    "countGenos", signature = "x",
    function(x, genos, pheno = NULL, level = NULL)
        standardGeneric("countGenos")
)

# variantsInSamples ----

setGeneric(
    "variantsInSamples", signature = c("vcf"),
    function(vcf, samples = 1:ncol(vcf), unique = FALSE)
        standardGeneric("variantsInSamples")
)

# vepInPhenoLevel ----

setGeneric(
    "vepInPhenoLevel", signature = c("vcf"),
    function(
        vcf, phenoCol, level, vepCol,
        unique = FALSE, facet = NULL)
        standardGeneric("vepInPhenoLevel")
)

# TVTBparam class ----

setGeneric(
    "TVTBparam", signature = "genos",
    function(
        genos, ranges = GRangesList(),
        aaf = "AAF", maf = "MAF", vep = "CSQ", bp = SerialParam(),
        svp = ScanVcfParam(which = reduce(unlist(ranges))))
        standardGeneric("TVTBparam")
)

setGeneric(
    "genos", signature = "x",
    function(x)
        standardGeneric("genos")
)

setGeneric(
    "genos<-", signature = c("x", "value"),
    function(x, value)
        standardGeneric("genos<-")
)

setGeneric(
    "aaf", signature = "x",
    function(x)
        standardGeneric("aaf")
)

setGeneric(
    "aaf<-", signature = c("x", "value"),
    function(x, value)
        standardGeneric("aaf<-")
)

setGeneric(
    "maf", signature = "x",
    function(x)
        standardGeneric("maf")
)

setGeneric(
    "maf<-", signature = c("x", "value"),
    function(x, value)
        standardGeneric("maf<-")
)

setGeneric(
    "het", signature = "x",
    function(x)
        standardGeneric("het")
)

setGeneric(
    "het<-", signature = c("x", "value"),
    function(x, value)
        standardGeneric("het<-")
)

setGeneric(
    "carrier", signature = "x",
    function(x)
        standardGeneric("carrier")
)

setGeneric(
    "vep", signature = "x",
    function(x)
        standardGeneric("vep")
)

setGeneric(
    "vep<-", signature = c("x", "value"),
    function(x, value)
        standardGeneric("vep<-")
)

setGeneric(
    "bp", signature = "x",
    function(x)
        standardGeneric("bp")
)

setGeneric(
    "bp<-", signature = c("x", "value"),
    function(x, value)
        standardGeneric("bp<-")
)

setGeneric(
    "suffix", signature = "x",
    function(x)
        standardGeneric("suffix")
)

setGeneric(
    "svp", signature = "x",
    function(x)
        standardGeneric("svp")
)

setGeneric(
    "svp<-", signature = c("x", "value"),
    function(x, value)
        standardGeneric("svp<-")
)
