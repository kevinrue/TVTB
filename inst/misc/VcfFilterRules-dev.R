

# Create a FilterRules object ---------------------------------------------

# library(S4Vectors)
#
# filts <- c("peaks", "promoters")
# filters <- FilterRules(filts)
# str(filters)
# elementMetadata(filters) <- list(a = "test")
#

# # define a VcfFilterRules class that extends FilterRules ----
# .valid.VcfFilterRules <- function(object){
#     return(TRUE)
# }
#
# VcfInfoRules <- setClass(
#     Class = "VcfInfoRules",
#
#     contains = "FilterRules",
#
#     validity = .valid.VcfFilterRules
# )
#
# VcfInfoRules # Class generator
# VcfInfoRules() # Empty object
# str(VcfInfoRules()) # Contains the same slots as FilterRules
#

# define the class generator ----
#
# setMethod(
#     f = "initialize",
#     signature = "VcfInfoRules",
#     definition = function(.Object, exprs = list(), ..., active = TRUE){
#
#         parentRules <- FilterRules(exprs = exprs, ..., active = TRUE)
#
#         .Object@listData <- slot(parentRules, "listData")
#         .Object@active <- slot(parentRules, "active")
#
#         validObject(.Object)
#         return(.Object)
#     }
# )

# import some variant in an ExpandedVCF ----

# VCF file
extdata <- file.path(system.file(package = "TVTB"), "extdata")
vcfFile <- file.path(extdata, "moderate.vcf")

# Parameters
library(TVTB)
tparam <- TVTBparam(
    genos = list(
        REF = "0|0",
        HET = c("0|1", "1|0"),
        ALT = "1|1"))

library(VariantAnnotation)
vcf <- readVcf(file = vcfFile)
vcf <- expand(x = vcf)
vcf <- addOverallFrequencies(vcf = vcf, param = tparam)

info(vcf)

# create an info filter ----
# 2-ways
# new(Class = "VcfInfoRules", exprs = list(INFOpass = expression(FILTER == "PASS")))
vif <- VcfInfoRules(exprs = list(
    common = expression(MAF > 0.1),
    alt = expression(ALT > 0)
    )) # Empty object
str(vif) # Contains the same slots as FilterRules

# Test info filter ----

eval(vif, vcf)
as.data.frame(evalSeparately(vif, vcf))

# Create a vep filter ----

vef <- VcfVepRules(exprs = list(
    missense = expression(Consequence %in% "missense_variant")
)) # Empty object
str(vef) # Contains the same slots as FilterRules

# Test vep filter ----

library(ensemblVEP)
csq <- parseCSQToGRanges(vcf, VCFRowID = rownames(vcf))
eval(expr = vef, envir = csq)
evalSeparately(expr = vef, envir = csq)

# Define an eval method for each type of filter ----

# setMethod(
#     f = "eval",
#     signature = c("VcfInfoRules", "ExpandedVCF"),
#     definition = function(expr, envir){
#         eval(expr = expr, envir = info(envir))
#     }
# )
#
# setMethod(
#     f = "eval",
#     signature = c(expr="VcfVepRules", envir="ExpandedVCF"),
#     definition = function(expr, envir){
#         # Extract VEP predictions
#         csq <- parseCSQToGRanges(x = envir, VCFRowID = rownames(envir))
#         # Does each prediction pass the filter
#         csqBoolean <- eval(expr = expr, envir = csq)
#         # Index of variants with 1+ prediction passing the filter
#         vcfPass <- mcols(csq)[csqBoolean,"VCFRowID", drop = TRUE]
#         # Does each variant have 1+ prediction passing the filter
#         vcfBoolean <- rep(FALSE, length(envir))
#         vcfBoolean[vcfPass] <- TRUE
#         return(vcfBoolean)
#     }
# )

# Test eval method for each type of filter ----

showMethods("eval")
eval(expr = vif, envir = vcf)
eval(expr = vef, envir = vcf)

# inherited
evalSeparately(expr = vif, envir = vcf)
evalSeparately(expr = vef, envir = vcf)

# Define a subset method for ExpandedVCF objects ----

# setMethod(
#     f = "subset",
#     signature = c(x="ExpandedVCF"),
#     definition = function(x, subset, ...){
#         # Apply filter
#         vcfBoolean <- eval(expr = subset, envir = x)
#         # Return records passing filter
#         return(x[vcfBoolean])
#     }
# )

# Test subset method for each type of filter ----

# subset(vcf, subset = vif)
# subset(vcf, subset = vef)

# Define a subsetByFilter method for each filter ----

# setMethod(
#     f = "subsetByFilter",
#     signature = c(x="ExpandedVCF", filter="VcfInfoRules"),
#     definition = function(x, filter, ...){
#         # Apply filter
#         vcfBoolean <- eval(expr = filter, envir = x)
#         # Return records passing filter
#         return(x[vcfBoolean])
#     }
# )
#
# setMethod(
#     f = "subsetByFilter",
#     signature = c(x="ExpandedVCF", filter="VcfVepRules"),
#     definition = function(x, filter, ...){
#         # Apply filter
#         vcfBoolean <- eval(expr = filter, envir = x)
#         # Return records passing filter
#         return(x[vcfBoolean])
#     }
# )

# Test subsetByFilter for each filter ----

subsetByFilter(x = vcf, filter = vif)
subsetByFilter(x = vcf, filter = vef)

# Extract the rules from a list of Vcf*Filters ----

# lf <- list(vif, vef)

# listLData <- unlist(lapply(lf, function(x){slot(x, "listData")}))
# lapply(listLData, class)
# vifLData <- slot(vif, "listData")
# unlist(vifLData)
# lapply(vif, class)
#
# unlist(lapply(lf, active))

# unlist(lapply(
#     X = lf,
#     FUN = function(x){
#             switch(
#                 class(x),
#                 VcfFixedRules = rep("fixed", length(x)),
#                 VcfInfoRules = rep("info", length(x)),
#                 VcfVepRules = rep("vep", length(x)),
#                 VcfFilterRules = slot(x, "type"),
#                 stop("Invalid filter")
#                 )
#         }
#     ))


# Combine Vcf*Rules in a VcfFilterRules
vfr <- VcfFilterRules(vif, vef)
str(vfr)
# Even VcfFilterRules can be combined into larger VcfFilterRules!
vfr.vef <- tryCatch(
    VcfFilterRules(vfr, vif),
    error = function(err){warning(err)})

# Test subsetByFilter for VcfFilterRules ----

vfr[1:2]
str(vfr[1:2])

eval(expr = vfr, envir = vcf)
evalSeparately(expr = vfr, envir = vcf)
subsetByFilter(x = vcf, filter = vfr)

# Type VcfFilterRules to fixed, info or vep (drop=TRUE) ----

VcfInfoRules(vfr[1:2]@listData)
str(VcfInfoRules(vfr[1:2]@listData))

vfr[1:2]
vfr[1:2, drop = TRUE]
tryCatch(
    vfr[2:3, drop = TRUE],
    error = function(err){warning(err)})
vfr[2:3, drop = FALSE]

vfr[[1]]
vfr[[1:2]]
vfr[[2:3]]
str(vfr[[2]])

str(vfr[[1]])
