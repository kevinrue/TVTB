context("dropInfo")

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
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)

vcfData <- vcf
extraInfo <- S4Vectors::DataFrame(FAKE = 1:length(vcfData))
info(vcfData) <- cbind(info(vcfData), extraInfo)

vcfHeader <- vcf
extraHeader <- S4Vectors::DataFrame(
    Number = 1,
    Type = "Integer",
    Description = "Completely fake",
    row.names = "FAKE"
)
info(header(vcfHeader)) <- rbind(info(header(vcfHeader)), extraHeader)

# Signatures ----

test_that("dropInfo() supports all signatures",{

    expect_s4_class(
        dropInfo(vcf = vcf),
        "VCF"
    )

})

# Drop a key present in both header and data ----

test_that("dropInfo() drops given key from header and data",{

    expect_s4_class(
        dropInfo(vcf = vcfHeader, key = "CSQ"),
        "ExpandedVCF"
    )

})

# Drop header ----

test_that("dropInfo() drops given key from header",{

    expect_s4_class(
        dropInfo(vcf = vcfHeader, key = "FAKE", slot = "header"),
        "ExpandedVCF"
    )

})

test_that("dropInfo() complains about key absent from header",{

    expect_message(
        dropInfo(vcf = vcfHeader, key = "FAK", slot = "header")
    )

})

# Drop data ----

test_that("dropInfo() drops given key from data",{

    expect_s4_class(
        dropInfo(vcf = vcfData, key = "FAKE", slot = "data"),
        "ExpandedVCF"
    )

})

test_that("dropInfo() complains about key absent from data",{

    expect_message(
        dropInfo(vcf = vcfData, key = "FAK", slot = "data")
    )

})
