# subsetVcf ----

setGeneric(
    name = "subsetVcf",
    signature = c("x", "filter", "param"),
    def = function(x, filter, param, ...)
        standardGeneric("subsetVcf")
)

# addOverallFrequencies ----

setGeneric(
    name = "addOverallFrequencies",
    signature = c("vcf", "param"),
    def = function(vcf, param, ...)
        standardGeneric("addOverallFrequencies")
)

# addPhenoLevelFrequencies ----

setGeneric(
    name = "addPhenoLevelFrequencies",
    signature = c("vcf", "param"),
    def = function(vcf, param, ...)
        standardGeneric("addPhenoLevelFrequencies")
)

# addPhenotypeFrequencies ----

setGeneric(
    name = "addFrequencies",
    signature = c("vcf", "phenos", "param"),
    def = function(vcf, phenos, param, ...)
        standardGeneric("addFrequencies")
)

# preprocessVariants ----

setGeneric(
    name = "preprocessVariants",
    signature = c("file", "param", "phenos"),
    def = function(
        file, param, phenos, ...)
        standardGeneric("preprocessVariants")
)

# addCountGenos ----

# Not used within the package anymore; left for user convenience
setGeneric(
    name = "addCountGenos",
    signature = c("vcf", "samples"),
    def = function(
        vcf, genos, key, description, samples = 1:ncol(vcf), force = FALSE)
        standardGeneric("addCountGenos")
)

# countGenos ----

setGeneric(
    name = "countGenos",
    signature = c("x","genos"),
    def = function(x, genos, ...)
        standardGeneric("countGenos")
)

# variantsInSamples ----

setGeneric(
    name = "variantsInSamples",
    signature = c("vcf", "samples", "param"),
    def = function(vcf, samples, param, ..., unique = FALSE)
        standardGeneric("variantsInSamples")
)

# csqInPhenoLevel ----

setGeneric(
    name = "csqInPhenoLevel",
    signature = c("vcf", "param"),
    def = function(
        vcf, phenoCol, level, csqCol, param, ...,
        unique = FALSE, facet = NULL)
        standardGeneric("csqInPhenoLevel")
)

# densityCsqByPhenotype ----

setGeneric(
    name = "densityCsqByPhenotype",
    signature = c("vcf", "param"),
    def = function(
        vcf, phenoCol, csqCol, param, ...,
        unique = FALSE, facet = NULL, plot = FALSE, popFreq = FALSE)
        standardGeneric("densityCsqByPhenotype")
)

# densityCsqInPhenoLevel ----

setGeneric(
    name = "densityCsqInPhenoLevel",
    signature = c("vcf", "param"),
    def = function(
        level, vcf, phenoCol, csqCol, param, ...,
        unique = FALSE, facet = NULL, plot = FALSE, popFreq = FALSE)
        standardGeneric("densityCsqInPhenoLevel")
)

# tabulateCsqByPhenotype ----

setGeneric(
    name = "tabulateCsqByPhenotype",
    signature = c("vcf", "param"),
    def = function(
        vcf, phenoCol, csqCol, param, ...,
        unique = FALSE, facet = NULL, plot = FALSE, percentage = FALSE)
        standardGeneric("tabulateCsqByPhenotype")
)

# tabulateCsqInPhenoLevel ----

setGeneric(
    name = "tabulateCsqInPhenoLevel",
    signature = c("vcf", "param"),
    def = function(
        level, vcf, phenoCol, csqCol, param, ...,
        unique = FALSE, facet = NULL, plot = FALSE, percentage = FALSE)
        standardGeneric("tabulateCsqInPhenoLevel")
)

# TVTBparam class ----

setGeneric(
    name = "TVTBparam",
    signature = "genos",
    def = function(
        ..., genos,
        ranges = GRanges(),
        aaf = "AAF", maf = "MAF", vep = "CSQ", bp = SerialParam())
        standardGeneric("TVTBparam")
)

setGeneric(
    name = "genos",
    signature = "x",
    def = function(x)
        standardGeneric("genos")
)

setGeneric(
    name = "genos<-",
    signature = c("x", "value"),
    def = function(x, value)
        standardGeneric("genos<-")
)

setGeneric(
    name = "aaf",
    signature = "x",
    def = function(x)
        standardGeneric("aaf")
)

setGeneric(
    name = "aaf<-",
    signature = c("x", "value"),
    def = function(x, value)
        standardGeneric("aaf<-")
)

setGeneric(
    name = "maf",
    signature = "x",
    def = function(x)
        standardGeneric("maf")
)

setGeneric(
    name = "maf<-",
    signature = c("x", "value"),
    def = function(x, value)
        standardGeneric("maf<-")
)

setGeneric(
    name = "hRef",
    signature = "x",
    def = function(x)
        standardGeneric("hRef")
)

setGeneric(
    name = "hRef<-",
    signature = c("x", "value"),
    def = function(x, value)
        standardGeneric("hRef<-")
)

setGeneric(
    name = "het",
    signature = "x",
    def = function(x)
        standardGeneric("het")
)

setGeneric(
    name = "het<-",
    signature = c("x", "value"),
    def = function(x, value)
        standardGeneric("het<-")
)

setGeneric(
    name = "hAlt",
    signature = "x",
    def = function(x)
        standardGeneric("hAlt")
)

setGeneric(
    name = "hAlt<-",
    signature = c("x", "value"),
    def = function(x, value)
        standardGeneric("hAlt<-")
)

setGeneric(
    name = "carrier",
    signature = "x",
    def = function(x)
        standardGeneric("carrier")
)

setGeneric(
    name = "carrier<-",
    signature = c("x", "value"),
    def = function(x, value)
        standardGeneric("carrier<-")
)

setGeneric(
    name = "vep",
    signature = "x",
    def = function(x)
        standardGeneric("vep")
)

setGeneric(
    name = "vep<-",
    signature = c("x", "value"),
    def = function(x, value)
        standardGeneric("vep<-")
)

setGeneric(
    name = "bp",
    signature = "x",
    def = function(x)
        standardGeneric("bp")
)

setGeneric(
    name = "bp<-",
    signature = c("x", "value"),
    def = function(x, value)
        standardGeneric("bp<-")
)

# VcfInfoFilter class ----

setGeneric(
    name = "VcfFixedFilter",
    signature = c("name", "condition", "value"),
    def = function(
        name, condition, value)
        standardGeneric("VcfFixedFilter")
)

setGeneric(
    name = "VcfInfoFilter",
    signature = c("name", "condition", "value"),
    def = function(
        name, condition, value)
        standardGeneric("VcfInfoFilter")
)

setGeneric(
    name = "VcfVepFilter",
    signature = c("name", "condition", "value"),
    def = function(
        name, condition, value)
        standardGeneric("VcfVepFilter")
)

setGeneric(
    name = "name",
    signature = "x",
    def = function(x)
        standardGeneric("name")
)

setGeneric(
    name = "name<-",
    signature = c("x", "value"),
    def = function(x, value)
        standardGeneric("name<-")
)

setGeneric(
    name = "condition",
    signature = "x",
    def = function(x)
        standardGeneric("condition")
)

setGeneric(
    name = "condition<-",
    signature = c("x", "value"),
    def = function(x, value)
        standardGeneric("condition<-")
)

setGeneric(
    name = "value",
    signature = "x",
    def = function(x)
        standardGeneric("value")
)

setGeneric(
    name = "value<-",
    signature = c("x", "value"),
    def = function(x, value)
        standardGeneric("value<-")
)

setGeneric(
    name = "filterType",
    signature = c("x"),
    def = function(x)
        standardGeneric("filterType")
)

# VcfFilterList class ----

setGeneric(
    name = "VcfFilterList",

    def = function(..., active = rep(TRUE, length(list(...))))
        standardGeneric("VcfFilterList")
)
