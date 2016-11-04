context("densityCsqByPhenotype")

# Settings ----

# VCF file
extdata <- file.path(system.file(package = "tSVE"), "extdata")
vcfFile <- file.path(extdata, "moderate.vcf")

# Good and bad phenotype files
phenoFile <- file.path(extdata, "moderate_pheno.txt")

tparam <- tSVEParam(
    genos = list(
        REF = "0|0",
        HET = c("0|1", "1|0"),
        ALT = "1|1"))

vcf <- preprocessVariants(
    file = vcfFile, param = tparam, phenos = phenoFile)

# Signatures ----

test_that("densityCsq* supports all signatures",{

    ## ByPhenotype()
    expect_is(
        densityCsqByPhenotype(
            vcf = vcf,
            phenoCol = "super_pop",
            csqCol = "CADD_PHRED",
            param = tparam),
        "data.frame"
    )

    expect_is(
        densityCsqByPhenotype(
            vcf = vcf,
            phenoCol = "super_pop",
            csqCol = "CADD_PHRED",
            alts = unlist(carrier(tparam))),
        "data.frame"
    )

    ## InPhenoLevel
    expect_is(
        densityCsqInPhenoLevel(
            level = "GBR",
            vcf = vcf,
            phenoCol = "pop",
            csqCol = "CADD_PHRED",
            param = tparam),
        "data.frame"
    )

    expect_is(
        densityCsqInPhenoLevel(
            level = "GBR",
            vcf = vcf,
            phenoCol = "pop",
            csqCol = "CADD_PHRED",
            alts = unlist(carrier(tparam))),
        "data.frame"
    )

})

# .checkAlts ----

test_that(".checkAlts catches invalid inputs", {

    expect_error(
        densityCsqByPhenotype(
            vcf = vcf,
            phenoCol = "super_pop",
            csqCol = "CADD_PHRED",
            alts = "0|1")
    )

})

# Arguments ----

test_that("plot & facet & popFreq argument work",{

    expect_s3_class(
        densityCsqByPhenotype(
            vcf = vcf,
            phenoCol = "super_pop",
            csqCol = "AMR_MAF",
            param = tparam,
            facet = "Feature", plot = TRUE, popFreq = TRUE),
        c("gg", "ggplot")
    )

    expect_s3_class(
        densityCsqInPhenoLevel(
            level = "AFR",
            vcf = vcf,
            phenoCol = "super_pop",
            csqCol = "AMR_MAF",
            param = tparam,
            facet = "Feature", plot = TRUE, popFreq = TRUE),
        c("gg", "ggplot")
    )

    # No variant in EUR population
    # expect_error(
    #     densityCsqInPhenoLevel(
    #         level = "EUR",
    #         vcf = vcf,
    #         phenoCol = "super_pop",
    #         csqCol = "AMR_MAF",
    #         param = tparam,
    #         facet = "Feature", plot = TRUE, popFreq = TRUE),
    #     c("gg", "ggplot")
    # )

})
