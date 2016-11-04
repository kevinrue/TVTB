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

fixedRules <- VcfFixedRules(
    exprs = list(
        pass = expression(FILTER == "PASS"),
        qual = function(env){env[,"QUAL"] > 20}
        ),
    active = c(TRUE, FALSE)
)

infoRules <- VcfInfoRules(
    exprs = list(
        common = expression(MAF > 0.1),
        alt = function(env){env[,"ALT"] > 0}
        ),
    active = c(TRUE, FALSE))

vepRules <- VcfVepRules(
    exprs = list(
        missense = expression(Consequence %in% c("missense_variant")),
        CADD = function(env){env[,"CADD_PHRED"] > 15}
        ),
    active = c(TRUE, FALSE))

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

# [[ ----

# NOTE: not critical as the parent method is called

test_that("[[ methods return valid values", {

    expect_that(fixedRules[[1]], is_a("expression"))
    expect_that(infoRules[[1]], is_a("expression"))
    expect_that(vepRules[[1]], is_a("expression"))
    expect_that(vcfRules[[1]], is_a("expression"))

    expect_s4_class(fixedRules[["qual"]], "FilterClosure")
    expect_s4_class(infoRules[["alt"]], "FilterClosure")
    expect_s4_class(vepRules[["CADD"]], "FilterClosure")
    expect_s4_class(vcfRules[["alt"]], "FilterClosure")

    expect_error(expect_s4_class(vcfRules[[1:2]]))

})

# [[<- ----

# NOTE: not critical as the parent method is called

allTrueFUN <- function(env){rep(TRUE, nrow(env))}

test_that("[[ methods perform valid replacement", {

    fixedRules[[1]] <- expression(PASS != "FAIL")
    fixedRules[["qual"]] <- allTrueFUN
    names(fixedRules)[2] <- "new_name"
    expect_that(fixedRules[[1]], is_a("expression"))
    expect_s4_class(fixedRules[["new_name"]], "FilterClosure")
    expect_equal(names(fixedRules)[2], "new_name")
    # replacing the expression of a filter reset active<-TRUE
    expect_equivalent(active(fixedRules), c(TRUE, TRUE))

    infoRules[[1]] <- expression(MAF > 0.05)
    infoRules[["alt"]] <- allTrueFUN
    names(infoRules)[2] <- "new_name"
    expect_that(infoRules[[1]], is_a("expression"))
    expect_s4_class(infoRules[["new_name"]], "FilterClosure")
    expect_equal(names(infoRules)[2], "new_name")
    expect_equivalent(active(infoRules), c(TRUE, TRUE))

    vepRules[[1]] <- expression(Consequence == c("missense_variant"))
    vepRules[["CADD"]] <- allTrueFUN
    names(vepRules)[2] <- "new_name"
    expect_that(vepRules[[1]], is_a("expression"))
    expect_s4_class(vepRules[["new_name"]], "FilterClosure")
    expect_equal(names(vepRules)[2], "new_name")
    expect_equivalent(active(vepRules), c(TRUE, TRUE))

    vcfRules[[1]] <- expression(QUAL %in% c("PASS"))
    vcfRules[["common"]] <- expression(MAF > 0.01)
    names(vcfRules)[4] <- "new_name"
    expect_that(vcfRules[[1]], is_a("expression"))
    expect_that(vcfRules[["common"]], is_a("expression"))
    expect_equal(names(vcfRules)[4], "new_name")
    # An active expression was replaced, therefore no active state changed
    expect_equivalent(active(vcfRules), rep(c(TRUE, FALSE), 3))
    expect_equal(type(vcfRules), rep(c("fixed", "info", "vep"), each = 2))

    expect_error(
        expect_s4_class(fixedRules[[1:2]] <- c(
            expression(A == 2),
            allTrueFUN))
    )

    expect_error(
        expect_s4_class(infoRules[[1:2]] <- c(
            expression(A == 2),
            allTrueFUN))
    )

    expect_error(
        expect_s4_class(vepRules[[1:2]] <- c(
            expression(A == 2),
            allTrueFUN))
    )

    expect_error(
        expect_s4_class(vcfRules[[1:2]] <- c(
            expression(A == 2),
            allTrueFUN))
    )

})

# [ ----

test_that("[ methods return valid values", {

    expect_equal(length(vcfRules[1:4]), 4)

    # Subset VcfFixedRules
    expect_s4_class(fixedRules[1:2], "VcfFixedRules")
    fixedSubset <- fixedRules["qual"]
    expect_s4_class(fixedSubset, "VcfFixedRules")
    expect_equivalent(active(fixedSubset), FALSE)

    # Subset VcfInfoRules
    expect_s4_class(infoRules[1:2], "VcfInfoRules")
    infoSubset <- infoRules["alt"]
    expect_s4_class(infoSubset, "VcfInfoRules")
    expect_equivalent(active(infoSubset), FALSE)

    # Subset VcfVepRules
    expect_s4_class(vepRules[1:2], "VcfVepRules")
    vep(vepRules) <- "ANN" # check that vep slot is transferred
    vepSubset <- vepRules["CADD"]
    expect_s4_class(vepSubset, "VcfVepRules")
    expect_equivalent(active(vepSubset), FALSE)
    expect_equal(vep(vepSubset), "ANN")

    # Subset VcfFilterRules
    ## subsetted objects are re-typed if possible
    ### Fixed
    expect_s4_class(vcfRules[1:2], "VcfFixedRules")
    expect_s4_class(vcfRules["pass"], "VcfFixedRules")
    ### Info
    expect_s4_class(vcfRules[3:4], "VcfInfoRules")
    expect_s4_class(vcfRules["common"], "VcfInfoRules")
    ### Vep (check that vep slot is transferred)
    expect_s4_class(vcfRules[5:6], "VcfVepRules")
    vep(vcfRules) <- "ANN2"
    vcfVepSubset <- vcfRules["CADD"]
    expect_s4_class(vcfVepSubset, "VcfVepRules")
    expect_equivalent(active(vcfVepSubset), FALSE)
    expect_equal(vep(vcfVepSubset), "ANN2")
    ### Vcf (check that vep & type slots are transferred)
    expect_s4_class(vcfRules[c(1,3,5)], "VcfFilterRules")
    vcfSubset <- vcfRules[c("CADD", "common")]
    expect_s4_class(vcfSubset, "VcfFilterRules")
    expect_equivalent(active(vcfSubset), c(FALSE, TRUE))
    expect_equal(vep(vcfSubset), "ANN2")

    # Subsetting is position-based, throw error if [row, column] given
    expect_error(
        expect_s4_class(vcfRules[1:2, 2])
    )

})

# [<- ----

test_that("[ methods perform valid replacement", {

    newFixedFilter <- VcfVepRules(exprs = list(
        filtSynonyms = expression(FILTER %in% c("PASS", "OK"))
    ))
    newInfoFilter <- VcfVepRules(exprs = list(
        altIsMinor = expression(AAF < 0.5)
    ))
    newVepFilter <- VcfVepRules(exprs = list(
        highImpact = expression(IMPACT == "HIGH")
    ))

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
