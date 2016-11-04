context("countGenos")

# Settings ----

# VCF file
vcfFile <- system.file("extdata", "moderate.vcf", package = "TVTB")

# Phenotype file
phenoFile <- system.file("extdata", "moderate_pheno.txt", package = "TVTB")
phenotypes <- S4Vectors::DataFrame(read.table(phenoFile, TRUE, row.names = 1))

# TVTB parameters
tparam <- TVTBparam(Genotypes("0|0", c("0|1", "1|0"), "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(
    vcfFile, param = tparam, colData = phenotypes)
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)

# countGenos() ---

test_that("countGenos() returns appropriate values",{

    expect_type(
        countGenos(vcf, het(tparam), "pop", "GBR"),
        "integer"
    )

    expect_type(
        countGenos(vcf, het(tparam)),
        "integer"
    )

})

# .checkPhenoLevel ----

test_that(".checkPhenoLevel() catches invalid inputs",{

    expect_error(countGenos(vcf, het(tparam), "missing", "GBR"))

    expect_error(countGenos(vcf, het(tparam), "pop", "missing"))

    expect_error(countGenos(vcf, het(tparam), level = "missing"))

})
