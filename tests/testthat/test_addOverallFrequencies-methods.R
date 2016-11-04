context("addOverallFrequencies")

# Settings ----

# VCF file
extdata <- file.path(system.file(package = "TVTB"), "extdata")
vcfFile <- file.path(extdata, "moderate.vcf")
phenoFile <- file.path(extdata, "moderate_pheno.txt")

tparam <- TVTBparam(
    genos = list(
        REF = "0|0",
        HET = c("0|1", "1|0"),
        ALT = "1|1"))

vcf <- preprocessVariants(
    file = vcfFile, param = tparam, phenos = phenoFile)

# Signatures ----

test_that("addOverallFrequencies supports all signatures",{

    # \alias{addOverallFrequencies,ExpandedVCF,TVTBparam-method}
    expect_s4_class(
        addOverallFrequencies(vcf = vcf, param = tparam),
        "ExpandedVCF"
    )

    # \alias{addOverallFrequencies,ExpandedVCF,missing-method}
    expect_s4_class(
        addOverallFrequencies(
            vcf = vcf,
            ref = unlist(hRef(tparam), use.names = FALSE),
            het = unlist(het(tparam), use.names = FALSE),
            alt = unlist(hAlt(tparam), use.names = FALSE)),
        "ExpandedVCF"
    )

})
