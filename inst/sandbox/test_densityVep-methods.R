context("densityVepByPhenotype")

# Settings ----

# VCF file
vcfFile <- system.file("extdata", "moderate.vcf", package = "TVTB")

# Phenotype file
phenoFile <- system.file("extdata", "moderate_pheno.txt", package = "TVTB")
phenotypes <- S4Vectors::DataFrame(read.table(phenoFile, TRUE, row.names = 1))

# TVTB parameters
tparam <- TVTBparam(Genotypes("0|0", c("0|1", "1|0"), "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(vcfFile)
colData(vcf) <- phenotypes
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)

# Signatures ----

test_that("densityVep* supports all signatures",{

    ## ByPhenotype()
    expect_is(
        densityVepByPhenotype(vcf, "super_pop", "CADD_PHRED", tparam),
        "data.frame"
    )

    ## InPhenoLevel
    expect_is(
        densityVepInPhenoLevel("GBR", vcf, "pop", "CADD_PHRED", tparam),
        "data.frame"
    )

})

# Arguments ----

test_that("plot & facet & pattern & layer arguments work",{

    expect_s3_class(
        densityVepByPhenotype(
            vcf, "super_pop", "AMR_MAF", tparam,
            facet = "Feature", plot = TRUE, pattern = ".*:(.*)"),
        c("gg", "ggplot")
    )

    expect_s3_class(
        densityVepInPhenoLevel(
            "AFR", vcf, "super_pop", "AMR_MAF", tparam,
            facet = "Feature", plot = TRUE, pattern = ".*:(.*)"),
        c("gg", "ggplot")
    )

    # No variant in EUR population
    # expect_error(
    #     densityVepInPhenoLevel(
    #         "EUR", vcf, "super_pop", "AMR_MAF", tparam,
    #         facet = "Feature", plot = TRUE, pattern = ".*:(.*)"),
    #     c("gg", "ggplot")
    # )

})
