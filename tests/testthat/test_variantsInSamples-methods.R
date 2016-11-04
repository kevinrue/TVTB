context("variantsInSamples")

# Settings ----

# VCF file
extdata <- file.path(system.file(package = "TVTB"), "extdata")
vcfFile <- file.path(extdata, "moderate.vcf")

# Phenotype file
phenoFile <- file.path(extdata, "moderate_pheno.txt")
phenotypes <- S4Vectors::DataFrame(
    read.table(file = phenoFile, header = TRUE, row.names = 1))

tparam <- TVTBparam(
    genos = list(
        REF = "0|0",
        HET = c("0|1", "1|0"),
        ALT = "1|1"))

vcf <- VariantAnnotation::readVcf(file = vcfFile)
colData(vcf) <- phenotypes
vcf <- VariantAnnotation::expand(vcf)

sampleIdx <- 1:(ncol(vcf))

# Signatures ----

test_that("variantsInSamples() supports all signatures",{

    ## ExpandedVCF,numeric,TVTBparam
    expect_type(
        variantsInSamples(
            vcf = vcf,
            samples = sampleIdx,
            param = tparam,
            unique = FALSE),
        "integer"
    )

    # # ExpandedVCF,numeric,missing / tested by higher functions TODO: change
    expect_type(
        variantsInSamples(
            vcf = vcf,
            samples = sampleIdx,
            alts = c("0|1", "1|0", "1|1"),
            unique = FALSE),
        "integer"
    )

    # ExpandedVCF,character,TVTBparam
    expect_type(
        variantsInSamples(
            vcf = vcf,
            samples = colnames(vcf)[sampleIdx],
            param = tparam,
            unique = FALSE),
        "integer"
    )

    # ExpandedVCF,character,missing
    expect_type(
        variantsInSamples(
            vcf = vcf,
            samples = colnames(vcf)[sampleIdx],
            alts = c("0|1", "1|0", "1|1"),
            unique = FALSE),
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
    #         alts = c("0|1", "1|0", "1|1"),
    #         unique = FALSE),
    #     "integer"
    # )

    # TRUE
    expect_type(
        variantsInSamples(
            vcf = vcf,
            samples = sampleIdx,
            alts = c("0|1", "1|0", "1|1"),
            unique = TRUE),
        "integer"
    )

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
