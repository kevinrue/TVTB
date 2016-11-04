context("addPhenoLevelFrequencies")

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

# Create a VCF object with a pre-existing INFO key

vcfInfoExist <- vcf

newInfoHeader <- DataFrame(
    Number = rep(1, 2),
    Type = "Integer",
    Description = "Pre-existing INFO field",
    row.names = c("MAF", "pop_GBR_MAF"))

newInfoData <- DataFrame(
    MAF = seq_along(vcfInfoExist),
    pop_GBR_MAF = rev(seq_along(vcfInfoExist))
)

info(header(vcfInfoExist)) <- rbind(info(header(vcfInfoExist)), newInfoHeader)
info(vcfInfoExist) <- cbind(info(vcfInfoExist), newInfoData)

# Signatures ----

test_that("addPhenoLevelFrequencies supports all signatures",{

    # \alias{addPhenoLevelFrequencies,ExpandedVCF,TVTBparam-method}
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

# .checkInputsPLF ----

test_that(".checkInputsPLF catches invalid inputs", {

    expect_error(
        addPhenoLevelFrequencies(
            vcf = vcf, pheno = "missing", level = "GBR", param = tparam)
    )

    expect_error(
        addPhenoLevelFrequencies(
            vcf = vcf, pheno = "pop", level = "missing", param = tparam)
    )

    expect_error(
        addPhenoLevelFrequencies(
            vcf = vcfInfoExist, pheno = "pop", level = "GBR", param = tparam)
    )

})
