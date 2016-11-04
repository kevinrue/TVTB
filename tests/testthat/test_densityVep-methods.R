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
colData(vcf) <- phenotypes
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)

# Signatures ----

test_that("densityVep* supports all signatures",{

    ## ByPhenotype()
    expect_is(
        densityVepByPhenotype(
            vcf = vcf,
            phenoCol = "super_pop",
            vepCol = "CADD_PHRED",
            param = tparam),
        "data.frame"
    )

    expect_is(
        densityVepByPhenotype(
            vcf = vcf,
            phenoCol = "super_pop",
            vepCol = "CADD_PHRED",
            alts = unlist(carrier(tparam))),
        "data.frame"
    )

    ## InPhenoLevel
    expect_is(
        densityVepInPhenoLevel(
            level = "GBR",
            vcf = vcf,
            phenoCol = "pop",
            vepCol = "CADD_PHRED",
            param = tparam),
        "data.frame"
    )

    expect_is(
        densityVepInPhenoLevel(
            level = "GBR",
            vcf = vcf,
            phenoCol = "pop",
            vepCol = "CADD_PHRED",
            alts = unlist(carrier(tparam))),
        "data.frame"
    )

})

# .checkAlts ----

test_that(".checkAlts catches invalid inputs", {

    expect_error(
        densityVepByPhenotype(
            vcf = vcf,
            phenoCol = "super_pop",
            vepCol = "CADD_PHRED",
            alts = "0|1")
    )

})

# Arguments ----

test_that("plot & facet & popFreq & layer arguments work",{

    expect_s3_class(
        densityVepByPhenotype(
            vcf = vcf,
            phenoCol = "super_pop",
            vepCol = "AMR_MAF",
            param = tparam,
            facet = "Feature", plot = TRUE, popFreq = TRUE),
        c("gg", "ggplot")
    )

    expect_s3_class(
        densityVepInPhenoLevel(
            level = "AFR",
            vcf = vcf,
            phenoCol = "super_pop",
            vepCol = "AMR_MAF",
            param = tparam,
            facet = "Feature", plot = TRUE, popFreq = TRUE),
        c("gg", "ggplot")
    )

    # No variant in EUR population
    # expect_error(
    #     densityVepInPhenoLevel(
    #         level = "EUR",
    #         vcf = vcf,
    #         phenoCol = "super_pop",
    #         vepCol = "AMR_MAF",
    #         param = tparam,
    #         facet = "Feature", plot = TRUE, popFreq = TRUE),
    #     c("gg", "ggplot")
    # )

})
