context("addFrequencies")

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
vcf <- addFrequencies(vcf, "super_pop")
# NOTE: Tamper with real values to avoid warning because all EUR_MAF = 0
# and cor() does not like invariant data
info(vcf)[,"super_pop_EUR_MAF"] <- runif(nrow(vcf), 0, 0.01)

# Signatures ----

test_that("all signatures work to completion", {

    expect_s3_class(
        pairsInfo(vcf, "MAF", "super_pop"),
        c("gg", "ggplot")
    )

})
