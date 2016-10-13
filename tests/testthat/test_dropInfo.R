context("dropInfo")

# Settings ----

# VCF file
vcfFile <- system.file("extdata", "moderate.vcf", package = "TVTB")

# TVTB parameters
tparam <- TVTBparam(Genotypes("0|0", c("0|1", "1|0"), "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(vcfFile, param = tparam)
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)

vcfHeader <- vcf
extraHeader <- S4Vectors::DataFrame(
    Number = 1, Type = "Integer", Description = "Completely fake",
    row.names = "FAKE"
)
VariantAnnotation::info(VariantAnnotation::header(vcfHeader)) <- rbind(
    VariantAnnotation::info(VariantAnnotation::header(vcfHeader)),
    extraHeader)

vcfData <- vcfHeader
extraInfo <- S4Vectors::DataFrame(FAKE = 1:length(vcfData))
VariantAnnotation::info(vcfData) <- cbind(
    VariantAnnotation::info(vcfData),
    extraInfo)


# Signatures ----

test_that("dropInfo() supports all signatures",{

    expect_s4_class(
        dropInfo(vcf),
        "VCF"
    )

})

# Drop a key present in both header and data ----

test_that("requested key are removed from header and data",{

    expect_s4_class(
        dropInfo(vcfHeader, "CSQ"),
        "ExpandedVCF"
    )

})

# Drop header ----

test_that("requested key are removed from header",{

    expect_s4_class(
        dropInfo(vcfHeader, "FAKE", "header"),
        "ExpandedVCF"
    )

})

test_that(".dropInfoHeader messages about absent keys",{

    expect_message(TVTB:::.dropInfoHeader(vcfHeader, "missing"))

})

# Drop data ----

test_that("requested key are removed from data",{

    expect_s4_class(
        dropInfo(vcfData, "FAKE", "data"),
        "ExpandedVCF"
    )

})

test_that(".dropInfoData messages about absent keys",{

    expect_message(dropInfo(vcfData, "missing", "data"))

})
