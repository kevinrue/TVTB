context("addOverallFrequencies")

# Settings ----

# VCF file
vcfFile <- system.file("extdata", "moderate.vcf", package = "TVTB")

# TVTB parameters
tparam <- TVTBparam(Genotypes("0|0", c("0|1", "1|0"), "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(
    vcfFile, param = tparam)
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)

# Create a VCF object with a pre-existing INFO key
vcfHeaderExist <- vcf
newInfoHeader <- S4Vectors::DataFrame(
    Number = rep(1, 2), Type = "Integer",
    Description = "Pre-existing INFO field",
    row.names = c("MAF", "pop_GBR_MAF"))
info(header(vcfHeaderExist)) <- rbind(
    info(header(vcfHeaderExist)),
    newInfoHeader)

# Create a VCF object with a pre-existing INFO key
vcfDataExist <- vcfHeaderExist
newInfoData <- S4Vectors::DataFrame(
    MAF = seq_along(vcfHeaderExist),
    pop_GBR_MAF = rev(seq_along(vcfHeaderExist))
)
info(vcfDataExist) <- cbind(info(vcfDataExist), newInfoData)

# Signatures ----

test_that("addOverallFrequencies supports all signatures",{

    # \alias{addOverallFrequencies,ExpandedVCF,TVTBparam-method}
    expect_s4_class(
        addOverallFrequencies(vcf),
        "ExpandedVCF"
    )

})

# Existing INFO header fields ----

test_that("Error thrown if INFO key cannot be overwritten", {

    expect_error(addOverallFrequencies(vcfHeaderExist))

    expect_error(addOverallFrequencies(vcfDataExist))

})

# Existing INFO data ----

test_that("Messages when overwriting INFO key", {

    expect_message(addOverallFrequencies(vcfHeaderExist, force = TRUE))

    expect_message(addOverallFrequencies(vcfDataExist, force = TRUE))

})
