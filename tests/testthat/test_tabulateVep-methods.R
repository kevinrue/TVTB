context("densityVepByPhenotype")

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

test_that("tabulateVep* supports all signatures",{

    # ExpandedVCF, TVTBparam
    expect_is(
        tabulateVepByPhenotype(
            vcf = vcf,
            phenoCol = "super_pop",
            vepCol = "Consequence",
            param = tparam,
            unique = TRUE),
        "data.frame"
    )

    expect_is(
        tabulateVepByPhenotype(
            vcf = vcf,
            phenoCol = "super_pop",
            vepCol = "Consequence",
            alts = c("0|1", "1|0", "1|1"),
            unique = TRUE),
        "data.frame"
    )

    ## Implicitely tested by *ByPhenotype
    expect_is(
        tabulateVepInPhenoLevel(
            level = "AFR",
            vcf = vcf,
            phenoCol = "super_pop",
            vepCol = "Consequence",
            param = tparam),
        "data.frame"
    )

    expect_is(
        tabulateVepInPhenoLevel(
            level = "AFR",
            vcf = vcf,
            phenoCol = "super_pop",
            vepCol = "Consequence",
            alts = c("0|1", "1|0", "1|1")),
        "data.frame"
    )

})

# Arguments ----

test_that("plot & facet & popFreq argument work",{

    expect_s3_class(
        tabulateVepByPhenotype(
            vcf = vcf, phenoCol = "super_pop",
            vepCol = "Consequence", param = tparam,
            facet = "Feature", plot = TRUE, percentage = TRUE),
        c("gg", "ggplot")
    )

    expect_s3_class(
        tabulateVepInPhenoLevel(
            level = "AFR", vcf = vcf,
            phenoCol = "super_pop",
            vepCol = "Consequence", param = tparam,
            facet = "Feature", plot = TRUE, percentage = TRUE),
        c("gg", "ggplot")
    )

    ## For 1% of extra coverage: plot=TRUE, percentage=FALSE,
    expect_s3_class(
        tabulateVepByPhenotype(
            vcf = vcf, phenoCol = "super_pop",
            vepCol = "Consequence", param = tparam,
            unique = FALSE, facet = "Feature",
            plot = TRUE),
        c("gg", "ggplot")
    )

    expect_s3_class(
        tabulateVepInPhenoLevel(
            level = "AFR", vcf = vcf,
            phenoCol = "super_pop",
            vepCol = "Consequence", param = tparam,
            unique = FALSE, facet = "Feature", plot = TRUE),
        c("gg", "ggplot")
    )
})

# test_that("errors are returned when ",{
#
#     expect_error(
#         tabulateVepInPhenoLevel(
#             level = "EUR", vcf = vcf,
#             phenoCol = "super_pop",
#             vepCol = "Consequence", param = tparam,
#             facet = "Feature", plot = TRUE, percentage = TRUE)
#     )
#
# })
