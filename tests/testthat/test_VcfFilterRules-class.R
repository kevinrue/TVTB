context("VcfFilterRules")

# Settings ----

# VCF file
vcfFile <- system.file("extdata", "moderate.vcf", package = "TVTB")

# TVTB parameters
tparam <- TVTBparam(Genotypes("0|0", c("0|1", "1|0"), "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(vcfFile, param = tparam)
vcf <- VariantAnnotation::expand(vcf)
vcf <- addOverallFrequencies(vcf, tparam)

# Test object ----

fixedRules <- VcfFixedRules(
    list(
        pass = expression(FILTER == "PASS"),
        qual = function(env){env[,"QUAL"] > 20}
        ),
    active = c(TRUE, FALSE)
)

infoRules <- VcfInfoRules(
    list(
        common = expression(MAF > 0.1),
        alt = function(env){env[,"ALT"] > 0}
        ),
    active = c(TRUE, FALSE))

vepRules <- VcfVepRules(
    list(
        missense = expression(Consequence %in% c("missense_variant")),
        CADD = function(env){env[,"CADD_PHRED"] > 15}
        ),
    active = c(TRUE, FALSE))

filterRules <- FilterRules(
    list(
        PASS = function(x) fixed(x)$FILTER == "PASS",
        QUAL = function(x) fixed(x)$QUAL > 20
    ),
    active = c(TRUE, FALSE)
)

vcfRules <- VcfFilterRules(fixedRules, infoRules, vepRules, filterRules)

filterNoVep <- VcfFilterRules(fixedRules, infoRules)

newFixedFilter <- VcfFixedRules(
    list(
        filtSynonyms = expression(FILTER %in% c("PASS", "OK")),
        fail = expression(FILTER=="FAIL")
        ),
    active = c(F, T))
newInfoFilter <- VcfInfoRules(
    list(
        altIsMinor = expression(AAF <= 0.5),
        altIsMajor = expression(AAF > 0.5)
        ),
    active = c(F, T))
newVepFilter <- VcfVepRules(
    list(
        highImpact = expression(IMPACT == "HIGH"),
        moderateImpact = expression(IMPACT == "MODERATE")
        ),
    active = c(F, T))
newVepSlot <- newVepFilter
vep(newVepSlot) <- "ANN"

# Constructors ----

test_that("Constructors produce a valid object",{

    expect_s4_class(fixedRules, "VcfFixedRules")

    expect_s4_class(infoRules, "VcfInfoRules")

    expect_s4_class(vepRules, "VcfVepRules")

    expect_s4_class(vcfRules, "VcfFilterRules")

    # Constructor including a VcfFilterRules object itself
    expect_s4_class(
        VcfFilterRules(VcfFilterRules(infoRules), vepRules),
        "VcfFilterRules"
    )

    # Incompatible vep slots
    vepANN <-  VcfVepRules(
        list(nonsense = expression(Fake == "NA")), vep = "ANN")
    expect_error(VcfFilterRules(vepRules, vepANN))

    # VcfFilterRules without vep rule
    vepNA <- VcfFilterRules(fixedRules, infoRules)
    expect_s4_class(
        vepNA,
        "VcfFilterRules"
    )
    expect_identical(vep(vepNA), NA_character_)

    expect_identical(filterNoVep@vep, NA_character_)
})

# Accessors ----

test_that("type accessor returns valid values", {

    type_VcfRules <- type(vcfRules)
    expect_type(type_VcfRules, "character")
    expect_length(type_VcfRules, length(vcfRules))
    expect_identical(names(type_VcfRules), names(type_VcfRules))

    type_VepRules <- type(vepRules)
    expect_identical(names(type_VepRules), names(vepRules))
    expect_equivalent(type_VepRules, rep("vep", length(vepRules)))

    type_InfoRules <- type(infoRules)
    expect_identical(names(type_InfoRules), names(type_InfoRules))
    expect_equivalent(type_InfoRules, rep("info", length(infoRules)))

    type_FixedRules <- type(fixedRules)
    expect_identical(names(type_FixedRules), names(type_FixedRules))
    expect_equivalent(type_FixedRules, rep("fixed", length(fixedRules)))

    type_FixedRules <- type(fixedRules)
    expect_identical(names(type_FixedRules), names(type_FixedRules))
    expect_equivalent(type_FixedRules, rep("fixed", length(fixedRules)))

})

test_that("vep accessor returns valid values", {

    vep_VcfRules <- vep(vcfRules)
    expect_type(vep_VcfRules, "character")
    expect_length(vep_VcfRules, 1)

    vep_VepRules <- vep(vepRules)
    expect_type(vep_VepRules, "character")
    expect_length(vep_VepRules, 1)

    vep_InfoRules <- vep(infoRules)
    expect_identical(vep_InfoRules, NA_character_)

    vep_FixedRules <- vep(fixedRules)
    expect_identical(vep_FixedRules, NA_character_)

    vep_FilterRules <- vep(filterRules)
    expect_identical(vep_FilterRules, NA_character_)

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
        eval(fixedRules, vcf),
        "logical")

    expect_type(
        eval(infoRules, vcf),
        "logical")

    expect_type(
        eval(vepRules, vcf),
        "logical")

    expect_type(
        eval(vcfRules, vcf),
        "logical")

})

test_that("evaluate empty filter rule is supported", {

    expect_identical(
        eval(VcfVepRules(), vcf),
        rep(TRUE, length(vcf))
    )

    expect_identical(
        eval(VcfFilterRules(), vcf),
        rep(TRUE, length(vcf))
    )

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

    expect_error(vcfRules[[1:2]])

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
    expect_equivalent(active(vcfRules), rep(c(TRUE, FALSE), 4))
    expect_equivalent(
        type(vcfRules),
        rep(c("fixed", "info", "vep", "filter"), each = 2)
    )

    expect_error(
        fixedRules[[1:2]] <- c(
            expression(A == 2),
            allTrueFUN)
    )

    expect_error(
        infoRules[[1:2]] <- c(
            expression(A == 2),
            allTrueFUN)
    )

    expect_error(
        vepRules[[1:2]] <- c(
            expression(A == 2),
            allTrueFUN)
    )

    expect_error(
        vcfRules[[1:2]] <- c(
            expression(A == 2),
            allTrueFUN)
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
    expect_error(fixedRules[1:2, 2])
    expect_error(infoRules[1:2, 2])
    expect_error(vepRules[1:2, 2])
    expect_error(vcfRules[1:2, 2])

    # Throw an error if name not found
    expect_error(vcfRules["missingFilter"])

})

# [<- ----

test_that("[<- can replace VcfFixedRules by index",{

    fixedRules[1:2] <- newFixedFilter

    expect_s4_class(fixedRules, "VcfFixedRules")
    expect_identical(fixedRules[[1]], newFixedFilter[[1]])
    expect_identical(fixedRules[[2]], newFixedFilter[[2]])

})

test_that("[<- can replace VcfFixedRules by name",{

    fixedRules["pass"] <- newFixedFilter["fail"]

    expect_s4_class(fixedRules, "VcfFixedRules")
    expect_identical(fixedRules[["pass"]], newFixedFilter[["fail"]])

})

test_that("[<- can replace VcfInfoRules by index",{

    infoRules[1:2] <- newInfoFilter

    expect_s4_class(infoRules, "VcfInfoRules")
    expect_identical(infoRules[[1]], newInfoFilter[[1]])
    expect_identical(infoRules[[2]], newInfoFilter[[2]])

})

test_that("[<- can replace VcfInfoRules by name",{

    infoRules["alt"] <- newInfoFilter["altIsMinor"]

    expect_s4_class(infoRules, "VcfInfoRules")
    expect_identical(infoRules[["alt"]], newInfoFilter[["altIsMinor"]])

})

test_that("[<- can replace VcfVepRules by index",{

    vepRules[1:2] <- newVepFilter

    expect_s4_class(vepRules, "VcfVepRules")
    expect_identical(vepRules[[1]], newVepFilter[[1]])
    expect_identical(vepRules[[2]], newVepFilter[[2]])

})

test_that("[<- can replace VcfVepRules by name",{

    vepRules["CADD"] <- newVepFilter["moderateImpact"]

    expect_s4_class(vepRules, "VcfVepRules")
    expect_identical(vepRules[["CADD"]], newVepFilter[["moderateImpact"]])

})

test_that("[<- can replace VcfFilterRules by index",{

    # Replacement is VcfFixedRules
    vcfRules[1:2] <- newFixedFilter

    expect_s4_class(vcfRules, "VcfFilterRules")
    expect_identical(vcfRules[[1]], newFixedFilter[[1]])
    expect_identical(vcfRules[[2]], newFixedFilter[[2]])

    # Replacement is VcfInfoRules
    vcfRules[3:4] <- newInfoFilter

    expect_s4_class(vcfRules, "VcfFilterRules")
    expect_identical(vcfRules[[3]], newInfoFilter[[1]])
    expect_identical(vcfRules[[4]], newInfoFilter[[2]])

    # Replacement is VcfVepRules
    vcfRules[5:6] <- newVepFilter

    expect_s4_class(vcfRules, "VcfFilterRules")
    expect_identical(vcfRules[[5]], newVepFilter[[1]])
    expect_identical(vcfRules[[6]], newVepFilter[[2]])

    # Replacement is VcfFilterRules
    vcfRules[5:6] <- VcfFilterRules(newVepFilter)

    expect_s4_class(vcfRules, "VcfFilterRules")
    expect_identical(vcfRules[[5]], newVepFilter[[1]])
    expect_identical(vcfRules[[6]], newVepFilter[[2]])
})

# Mimic behaviour of base::list
# TODO: requires fix in S4Vector::Vector
test_that("[<-NULL removes rules", {

    # fixedRules[1] <- NULL
    # expect_equal(length(fixedRules), 1) # 2 - 1

    # infoRules[1] <- NULL
    # expect_equal(length(fixedRules), 1) # 2 - 1

    # vepRules[1] <- NULL
    # expect_equal(length(fixedRules), 1) # 2 - 1

    # vcfRules[1] <- NULL
    # expect_equal(length(vcfRules), 5) # 6 - 1

    # vcfRules[1:3] <- NULL
    # expect_s4_class(vcfRules, "VcfFilterRules")
})

test_that("[<- requires compatible vep slots", {

    # x=VcfVepRules
    expect_error(vepRules[1:2] <- newVepSlot)

    # x=VcfFilterRules
    expect_error(vcfRules[1:2] <- newVepSlot)

})


# TODO: requires fix in S4Vector::FilterRules to subset active slot
test_that("[ requires valid length of values", {

    # currently produces error & warning
    # expect_warning(fixedRules[1] <- newFixedFilter) # 1 vs 2
    # expect_warning(infoRules[1] <- newInfoFilter) # 1 vs 2
    # expect_warning(vepRules[1] <- newVepFilter) # 1 vs 2
    # expect_warning(vcfRules[1:3] <- newFixedFilter) # 3 vs 2

})

test_that("[ requires value of a valid class", {

    expect_error(fixedRules[1] <- newFixedFilter[[1]])
    expect_error(infoRules[1] <- newInfoFilter[[1]])
    expect_error(vepRules[1] <- newVepFilter[[1]])
    expect_error(fixedRules[1] <- newVepFilter[[1]])

})

test_that("[ throws an error if columns indices are given", {

    expect_error(fixedRules[1:2, 2] <- NULL)
    expect_error(infoRules[1:2, 2] <- NULL)
    expect_error(vepRules[1:2, 2] <- NULL)
    expect_error(vcfRules[1:2, 2] <- NULL)

})

# coerce method ----

test_that("as method coerces objects", {

    expect_s4_class(
        as(fixedRules, "VcfFilterRules"),
        "VcfFilterRules"
    )

    expect_s4_class(
        as(infoRules, "VcfFilterRules"),
        "VcfFilterRules"
    )

    expect_s4_class(
        as(vepRules, "VcfFilterRules"),
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

# append method ----

test_that("append method combine values", {

    expect_equal(
        length(append(fixedRules, infoRules)),
        length(fixedRules) + length(infoRules)
    )

    expect_equal(
        length(append(infoRules, vepRules)),
        length(infoRules) + length(vepRules)
    )

    expect_equal(
        length(append(vepRules, fixedRules)),
        length(vepRules) + length(fixedRules)
    )

    expect_equal(
        length(append(filterNoVep, vepRules)),
        length(filterNoVep) + length(vepRules)
    )

})

# .dropVcfFilterRules method ----

test_that("VcfFilterRules of length 0 cannot be down-typed", {

    expect_identical(
        TVTB:::.dropVcfFilterRules(TVTB::VcfFilterRules()),
        TVTB::VcfFilterRules()
    )

})
