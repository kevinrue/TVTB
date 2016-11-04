context("variantsInSamples")

# Settings ----

# VCF file
vcfFile <- system.file("extdata", "moderate.vcf", package = "TVTB")

# TVTB parameters
tparam <- TVTBparam(Genotypes("0|0", c("0|1", "1|0"), "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(vcfFile, param = tparam)
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)

sampleIdx <- 1:(ncol(vcf))

# Signatures ----

test_that("all signatures are supported",{

    ## samples: numeric
    expect_type(
        variantsInSamples(vcf, sampleIdx, unique = FALSE),
        "integer"
    )

})

test_that("unique=TRUE is supported", {

    expect_type(
        variantsInSamples(
            vcf, colnames(vcf)[sampleIdx], unique = TRUE),
        "integer"
    )

})
