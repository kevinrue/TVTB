context("variantsInSamples")

# Settings ----

# VCF file
vcfFile <- system.file("extdata", "moderate.vcf", package = "TVTB")

# Import variants
vcf <- VariantAnnotation::readVcf(vcfFile)

# Signatures ----

test_that("new TVTBparam added in metadata if missing", {

    # Add
    expect_message(
        "TVTBparam" %in% names(metadata(autodetectGenotypes(vcf))),
        "created"
    )

})
