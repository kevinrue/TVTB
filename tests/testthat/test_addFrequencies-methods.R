context("addFrequencies")

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

test_that("addFrequencies supports all signatures",{

    # \alias{addFrequencies,ExpandedVCF,list,tSVEParam-method}
    expect_s4_class(
        addFrequencies(
            vcf = vcf, phenos = list(pop = "GBR"), param = tparam),
        "ExpandedVCF"
    )
    # \alias{addFrequencies,ExpandedVCF,list,missing-method}
    expect_s4_class(
        addFrequencies(
            vcf = vcf, phenos = list(pop = "GBR"),
            ref = unlist(hRef(tparam), use.names = FALSE),
            het = unlist(het(tparam), use.names = FALSE),
            alt = unlist(hAlt(tparam), use.names = FALSE)),
        "ExpandedVCF"
    )

    # \alias{addFrequencies,ExpandedVCF,character,tSVEParam-method}
    expect_s4_class(
        addFrequencies(
            vcf = vcf, phenos = "gender", param = tparam),
        "ExpandedVCF"
    )
    # \alias{addFrequencies,ExpandedVCF,character,missing-method}
    expect_s4_class(
        addFrequencies(
            vcf = vcf, phenos = "gender",
            ref = unlist(hRef(tparam), use.names = FALSE),
            het = unlist(het(tparam), use.names = FALSE),
            alt = unlist(hAlt(tparam), use.names = FALSE)),
        "ExpandedVCF"
    )

    # \alias{addFrequencies,ExpandedVCF,missing,tSVEParam-method}
    expect_s4_class(
        addFrequencies(
            vcf = vcf, param = tparam),
        "ExpandedVCF"
    )
    # \alias{addFrequencies,ExpandedVCF,missing,missing-method}
    expect_s4_class(
        addFrequencies(
            vcf = vcf,
            ref = unlist(hRef(tparam), use.names = FALSE),
            het = unlist(het(tparam), use.names = FALSE),
            alt = unlist(hAlt(tparam), use.names = FALSE)),
        "ExpandedVCF"
    )
})
