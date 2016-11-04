context("variantsInSamples")

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

sampleIdx <- 1:(ncol(vcf))

# Signatures ----

test_that("variantsInSamples() supports all signatures",{

    ## samples: numeric
    expect_type(
        variantsInSamples(
            vcf = vcf,
            samples = sampleIdx,
            param = tparam,
            unique = FALSE),
        "integer"
    )

    # ExpandedVCF,TVTBparam
    expect_type(
        variantsInSamples(
            vcf = vcf,
            samples = colnames(vcf)[sampleIdx],
            param = tparam,
            unique = TRUE),
        "integer"
    )

})

# Argument: unique ----

test_that("unique=TRUE is supported",{

    # FALSE / tested above
    # expect_type(
    #     variantsInSamples(
    #         vcf = vcf,
    #         samples = sampleIdx,
    #         param = tparam,
    #         unique = FALSE),
    #     "integer"
    # )

    # TRUE / tested above
    # expect_type(
    #     variantsInSamples(
    #         vcf = vcf,
    #         samples = sampleIdx,
    #         param = tparam,
    #         unique = TRUE),
    #     "integer"
    # )

})

# Argument: unique ----

test_that("two or more alternate genotypes are required",{

    expect_error(
        variantsInSamples(
            vcf = vcf,
            samples = sampleIdx,
            alts = c("1|1"),
            unique = FALSE)
    )

    expect_error(
        variantsInSamples(
            vcf = vcf,
            samples = colnames(vcf)[sampleIdx],
            alts = c("1|1"),
            unique = FALSE)
    )

})
