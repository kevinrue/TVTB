context("addOverallFrequencies")

# Settings ----

# VCF file
extdata <- file.path(system.file(package = "TVTB"), "extdata")
vcfFile <- file.path(extdata, "moderate.vcf")

# TVTB parameters
tparam <- TVTBparam(
    genos = list(
        REF = "0|0",
        HET = c("0|1", "1|0"),
        ALT = "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(file = vcfFile)
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)

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
