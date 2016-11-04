context("vepInPhenoLevel")

# Settings ----

# VCF file
extdata <- file.path(system.file(package = "TVTB"), "extdata")
vcfFile <- file.path(extdata, "moderate.vcf")

# Phenotype file
phenoFile <- file.path(extdata, "moderate_pheno.txt")
phenotypes <- S4Vectors::DataFrame(
    read.table(file = phenoFile, header = TRUE, row.names = 1))

# TVTB parameters
tparam <- TVTBparam(
    genos = list(
        REF = "0|0",
        HET = c("0|1", "1|0"),
        ALT = "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(file = vcfFile)
colData(vcf) <- phenotypes
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)

# Signatures ----

test_that("vepInPhenoLevel() supports all signatures",{

    # ExpandedVCF, / implicitely tested by higher functions TODO:
    expect_is(
        vepInPhenoLevel(
            vcf = vcf, phenoCol = "super_pop", level = "EUR",
            vepCol = "CADD_PHRED", param = tparam,
            unique = FALSE, facet = NULL),
        "data.frame"
    )

})
