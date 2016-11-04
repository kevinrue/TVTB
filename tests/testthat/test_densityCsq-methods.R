context("densityCsqByPhenotype")

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
