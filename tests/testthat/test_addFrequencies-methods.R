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

# Signatures ----

test_that("addFrequencies supports all signatures",{

    # \alias{addFrequencies,ExpandedVCF,list-method}
    expect_s4_class(
        addFrequencies(vcf, list(pop = "GBR")),
        "ExpandedVCF"
    )

    # \alias{addFrequencies,ExpandedVCF,character-method}
    expect_s4_class(
        addFrequencies(vcf, "gender"),
        "ExpandedVCF"
    )

    # \alias{addFrequencies,ExpandedVCF,missing-method}
    expect_s4_class(
        addFrequencies(vcf),
        "ExpandedVCF"
    )

})
