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
# Add phenotype information necessary for the demo
colData(vcf) <- phenotypes
# Separate multi-allelic records into bi-allelic records
vcf <- VariantAnnotation::expand(x = vcf, row.names = TRUE)
# Disambiguate row.names from multi-allelic records
rownames(vcf) <- paste(rownames(vcf), mcols(vcf)[,"ALT"], sep = "_")

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

    # ExpandedVCF, character
    expect_is(
        vepInPhenoLevel(
            vcf = vcf, phenoCol = "super_pop", level = "EUR",
            vepCol = "CADD_PHRED",
            alts = unlist(carrier(tparam)),
            unique = FALSE, facet = NULL),
        "data.frame"
    )

})

# Invalid input ----

test_that("vepInPhenoLevel() catches invalid input",{

    expect_error(
        vepInPhenoLevel(
            vcf = vcf, phenoCol = "super_pop", level = "EUR",
            vepCol = "CADD_PHRED", param = tparam,
            het = NULL, alt = "1|1",
            unique = FALSE, facet = NULL)
    )

    # ExpandedVCF, character
    expect_error(
        vepInPhenoLevel(
            vcf = vcf, phenoCol = "super_pop", level = "EUR",
            vepCol = "CADD_PHRED",
            alts = "0|1",
            unique = FALSE, facet = NULL)
    )

})
