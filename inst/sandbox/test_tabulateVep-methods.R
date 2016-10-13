context("densityVepByPhenotype")

# Settings ----

# VCF file
vcfFile <- system.file("extdata", "moderate.vcf", package = "TVTB")

# Phenotype file
phenoFile <- system.file("extdata", "moderate_pheno.txt", package = "TVTB")
phenotypes <- S4Vectors::DataFrame(
    read.table(phenoFile, header = TRUE, row.names = 1))

# TVTB parameters
tparam <- TVTBparam(Genotypes("0|0", c("0|1", "1|0"), "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(vcfFile)
colData(vcf) <- phenotypes
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)

# Signatures ----

test_that("tabulateVep* supports all signatures",{

    # ExpandedVCF, TVTBparam
    expect_is(
        tabulateVepByPhenotype(
            vcf, "super_pop", "Consequence", tparam,
            unique = TRUE, facet = "Feature"),
        "data.frame"
    )

    ## Implicitely tested by *ByPhenotype
    expect_is(
        tabulateVepInPhenoLevel(
            "AFR", vcf, "super_pop", "Consequence", tparam,
            facet = "Feature"),
        "data.frame"
    )

})

# Arguments ----

test_that("plot & facet & popFreq argument work",{

    expect_s3_class(
        tabulateVepByPhenotype(
            vcf,"super_pop", "Consequence", tparam,
            facet = "Feature", plot = TRUE, percentage = TRUE),
        c("gg", "ggplot")
    )

    expect_s3_class(
        tabulateVepInPhenoLevel(
            "AFR", vcf, "super_pop", "Consequence", tparam,
            facet = "Feature", plot = TRUE, percentage = TRUE),
        c("gg", "ggplot")
    )

    ## For 1% of extra coverage: plot=TRUE, percentage=FALSE,
    expect_s3_class(
        tabulateVepByPhenotype(
            vcf, "super_pop", "Consequence", tparam,
            unique = FALSE, facet = "Feature", plot = TRUE),
        c("gg", "ggplot")
    )

    expect_s3_class(
        tabulateVepInPhenoLevel(
            "AFR", vcf, "super_pop", "Consequence", tparam,
            unique = FALSE, facet = "Feature", plot = TRUE),
        c("gg", "ggplot")
    )

})

# test_that("errors are returned when ",{
#
#     expect_error(
#         tabulateVepInPhenoLevel(
#             "EUR", vcf, "super_pop", "Consequence", tparam,
#             facet = "Feature", plot = TRUE, percentage = TRUE)
#     )
#
# })
