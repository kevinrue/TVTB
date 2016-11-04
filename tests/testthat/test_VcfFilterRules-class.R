context("VcfFilterRules")

# Settings ----

# VCF file
extdata <- file.path(system.file(package = "TVTB"), "extdata")
vcfFile <- file.path(extdata, "moderate.vcf")

# TVTB parameters
tparam <- TVTBparam(
    genos = list(
        REF = "0|0",
        HET = c("0|1", "1|0"),
        ALT = "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(file = vcfFile)
vcf <- VariantAnnotation::expand(vcf)
vcf <- addOverallFrequencies(vcf = vcf, param = tparam)

# Test object ----

fixedRules <- VcfFixedRules(exprs = list(
    pass = expression(FILTER == "PASS"),
    qual = expression(QUAL > 20)
))

infoRules <- VcfInfoRules(exprs = list(
    common = expression(MAF > 0.1), # minor allele frequency
    alt = expression(ALT > 0) # count of alternative homozygotes
))

vepRules <- VcfVepRules(exprs = list(
    missense = expression(Consequence %in% c("missense_variant")),
    CADD = expression(CADD_PHRED > 15)
))

vcfRules <- VcfFilterRules(fixedRules, infoRules, vepRules)
filterNoVep <- VcfFilterRules(fixedRules, infoRules)

# Constructors ----

test_that("Constructors produce a valid object",{

    expect_s4_class(fixedRules, "VcfFixedRules")

    expect_s4_class(infoRules, "VcfInfoRules")

    expect_s4_class(vepRules, "VcfVepRules")

    expect_s4_class(vcfRules, "VcfFilterRules")

    vepANN <-  VcfVepRules(exprs = list(
        nonsense = expression(Fake == "NA")), vep = "ANN")
    expect_error(VcfFilterRules(vepRules, vepANN))

    expect_identical(filterNoVep@vep, NA_character_)
})

# Accessors ----

test_that("Accessors return valid values", {

    expect_type(
        type(vcfRules),
        "character")

})

# Setters ----

# test_that("Setters return valid values", {
#
#     expect_type(
#         vep(vcfRules) <- "ANN",
#         "character")
#
# })

# eval method ----

test_that("eval method return valid values", {

    expect_type(
        eval(expr = fixedRules, envir = vcf),
        "logical")

    expect_type(
        eval(expr = infoRules, envir = vcf),
        "logical")

    expect_type(
        eval(expr = vepRules, envir = vcf),
        "logical")

    expect_type(
        eval(expr = vcfRules, envir = vcf),
        "logical")

})

# Subset VcfFilterRules ----

test_that("[ and [[ methods return valid values", {

    vcfRulesMix <- vcfRules[seq(from=2, to=6, by=2)]

    expect_equal(length(vcfRulesMix), 3)
    expect_s4_class(vcfRules[1:2], "VcfFixedRules")
    expect_s4_class(vcfRules[3:4], "VcfInfoRules")
    expect_s4_class(vcfRules[5:6], "VcfVepRules")
    expect_s4_class(vcfRules["CADD"], "VcfVepRules")

    vcfRulesSingle <- vcfRules[[4]]
    expect_equal(length(vcfRulesSingle), 1)
    expect_s4_class(vcfRulesSingle, "VcfInfoRules")

    expect_error(
        expect_s4_class(vcfRules[1:2, 2])
    )

    expect_error(
        expect_s4_class(vcfRules[[1:2, 2]])
    )

})

# Replace VcfFilterRules element ----

test_that("[ and [[ methods assign values", {

    newFixedFilter <- VcfVepRules(exprs = list(
        filtSynonyms = expression(FILTER %in% c("PASS", "OK"))
    ))
    newInfoFilter <- VcfVepRules(exprs = list(
        altIsMinor = expression(AAF < 0.5)
    ))
    newVepFilter <- VcfVepRules(exprs = list(
        highImpact = expression(IMPACT == "HIGH")
    ))

    fixedRules[[1]] <- newFixedFilter
    expect_identical(fixedRules[[1]], newFixedFilter)

    infoRules[[1]] <- newInfoFilter
    expect_identical(infoRules[[1]], newInfoFilter)

    vepRules[[1]] <- newVepFilter
    expect_identical(vepRules[[1]], newVepFilter)

    vcfRules[[1]] <- newVepFilter
    expect_identical(vcfRules[[1]], newVepFilter)

    vcfRules[["highImpact"]] <- newVepFilter
    expect_identical(vcfRules[[1]], newVepFilter)

    expect_error(
        expect_s4_class(fixedRules[[1, 2]] <- newFixedFilter)
    )

    expect_error(
        expect_s4_class(infoRules[[1, 2]] <- newInfoFilter)
    )

    expect_error(
        expect_s4_class(vepRules[[1, 2]] <- newVepFilter)
    )

    expect_error(
        expect_s4_class(vcfRules[[1, 2]] <- newVepFilter)
    )

})

# coerce method ----

test_that("as method coerces objects", {

    expect_s4_class(
        as(object = fixedRules, Class = "VcfFilterRules"),
        "VcfFilterRules"
    )

    expect_s4_class(
        as(object = infoRules, Class = "VcfFilterRules"),
        "VcfFilterRules"
    )

    expect_s4_class(
        as(object = vepRules, Class = "VcfFilterRules"),
        "VcfFilterRules"
    )

})

# c method ----

test_that("c method combine values", {

    expect_equal(
        length(c(fixedRules, infoRules)),
        length(fixedRules) + length(infoRules)
    )

    expect_equal(
        length(c(infoRules, vepRules)),
        length(infoRules) + length(vepRules)
    )

    expect_equal(
        length(c(vepRules, fixedRules)),
        length(vepRules) + length(fixedRules)
    )

    expect_equal(
        length(c(filterNoVep, vepRules)),
        length(filterNoVep) + length(vepRules)
    )

})
