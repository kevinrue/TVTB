# dropInfo ----

setGeneric(
    name = "dropInfo",
    signature = c("vcf"),
    def = function(vcf, key = NULL, slot = "both")
        standardGeneric("dropInfo")
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
    def = function(vcf, pheno, level, param, ...)
        standardGeneric("addPhenoLevelFrequencies")
)

# addPhenotypeFrequencies ----

setGeneric(
    name = "addFrequencies",
    signature = c("vcf", "phenos", "param"),
    def = function(vcf, phenos, param, ...)
        standardGeneric("addFrequencies")
)

# addCountGenos ----

# Not used within the package anymore; left for user convenience
setGeneric(
    name = "addCountGenos",
    signature = c("vcf"),
    def = function(
        vcf, genos, key, description, samples = 1:ncol(vcf), force = FALSE)
        standardGeneric("addCountGenos")
)

# countGenos ----

setGeneric(
    name = "countGenos",
    signature = c("x"),
    def = function(x, genos, pheno = NULL, level = NULL)
        standardGeneric("countGenos")
)

# variantsInSamples ----

setGeneric(
    name = "variantsInSamples",
    signature = c("vcf", "samples", "param"),
    def = function(vcf, samples, param, ..., unique = FALSE)
        standardGeneric("variantsInSamples")
)

# vepInPhenoLevel ----

setGeneric(
    name = "vepInPhenoLevel",
    signature = c("vcf", "param"),
    def = function(
        vcf, phenoCol, level, vepCol, param, ...,
        unique = FALSE, facet = NULL)
        standardGeneric("vepInPhenoLevel")
)

# densityVepByPhenotype ----

setGeneric(
    name = "densityVepByPhenotype",
    signature = c("vcf", "param"),
    def = function(
        vcf, phenoCol, vepCol, param, ..., filter = VcfFilterRules(),
        unique = FALSE, facet = NULL, plot = FALSE, pattern = NULL,
        layer = "density+dotplot")
        standardGeneric("densityVepByPhenotype")
)

# densityVepInPhenoLevel ----

setGeneric(
    name = "densityVepInPhenoLevel",
    signature = c("vcf", "param"),
    def = function(
        level, vcf, phenoCol, vepCol, param, ..., filter = VcfFilterRules(),
        unique = FALSE, facet = NULL, plot = FALSE, pattern = NULL,
        layer = "density+dotplot")
        standardGeneric("densityVepInPhenoLevel")
)

# tabulateVepByPhenotype ----

setGeneric(
    name = "tabulateVepByPhenotype",
    signature = c("vcf", "param"),
    def = function(
        vcf, phenoCol, vepCol, param, ..., filter = VcfFilterRules(),
        unique = FALSE, facet = NULL, plot = FALSE, percentage = FALSE)
        standardGeneric("tabulateVepByPhenotype")
)

# tabulateVepInPhenoLevel ----

setGeneric(
    name = "tabulateVepInPhenoLevel",
    signature = c("vcf", "param"),
    def = function(
        level, vcf, phenoCol, vepCol, param, ..., filter = VcfFilterRules(),
        unique = FALSE, facet = NULL, plot = FALSE, percentage = FALSE)
        standardGeneric("tabulateVepInPhenoLevel")
)

# TVTBparam class ----

setGeneric(
    name = "TVTBparam",
    signature = "genos",
    def = function(
        genos,
        ranges = GRangesList(),
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
