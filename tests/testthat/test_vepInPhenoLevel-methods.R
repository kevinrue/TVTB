context("vepInPhenoLevel")

# Settings ----

# VCF file
vcfFile <- system.file("extdata", "moderate.vcf", package = "TVTB")

# Phenotype file
phenoFile <- system.file("extdata", "moderate_pheno.txt", package = "TVTB")
phenotypes <- S4Vectors::DataFrame(
    read.table(file = phenoFile, header = TRUE, row.names = 1))

# TVTB parameters
tparam <- TVTBparam(Genotypes("0|0", c("0|1", "1|0"), "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(
    vcfFile, param = tparam, colData = phenotypes)
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)

# Signatures ----

test_that("vepInPhenoLevel() supports all signatures",{

    # ExpandedVCF, / implicitely tested by higher functions TODO:
    expect_is(
        vepInPhenoLevel(
            vcf, "super_pop", "AFR", "CADD_PHRED",
            unique = FALSE, facet = NULL),
        "data.frame"
    )

})

# No data to return after filtering ----

test_that("vepInPhenoLevel() supports all signatures",{

    # ExpandedVCF, / implicitely tested by higher functions TODO:
    expect_is(
        vepInPhenoLevel(
            vcf, "super_pop", "EUR", "CADD_PHRED",
            unique = FALSE, facet = NULL),
        "data.frame"
    )

})

