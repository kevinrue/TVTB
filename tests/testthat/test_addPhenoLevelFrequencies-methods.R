context("addPhenoLevelFrequencies")

# Settings ----

# VCF file
extdata <- file.path(system.file(package = "tSVE"), "extdata")
vcfFile <- file.path(extdata, "moderate.vcf")
phenoFile <- file.path(extdata, "moderate_pheno.txt")

tparam <- tSVEParam(
    genos = list(
        REF = "0|0",
        HET = c("0|1", "1|0"),
        ALT = "1|1"))

vcf <- preprocessVariants(
    file = vcfFile, param = tparam, phenos = phenoFile)

# Signatures ----

test_that("addPhenoLevelFrequencies supports all signatures",{

    # \alias{addPhenoLevelFrequencies,ExpandedVCF,tSVEParam-method}
    expect_s4_class(
        addPhenoLevelFrequencies(
            vcf = vcf, pheno = "pop", level = "GBR", param = tparam),
        "ExpandedVCF"
    )

    # \alias{addPhenoLevelFrequencies,ExpandedVCF,missing-method}
    expect_s4_class(
        addPhenoLevelFrequencies(
            vcf = vcf, pheno = "pop", level = "GBR",
            ref = unlist(hRef(tparam), use.names = FALSE),
            het = unlist(het(tparam), use.names = FALSE),
            alt = unlist(hAlt(tparam), use.names = FALSE)),
        "ExpandedVCF"
    )

})
