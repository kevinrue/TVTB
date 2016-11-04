context("addFrequencies")

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

# Create a VCF object with a pre-existing INFO key

vcfHeaderExist <- vcf
newInfoHeader <- DataFrame(
    Number = rep(1, 2),
    Type = "Integer",
    Description = "Pre-existing INFO field",
    row.names = c("MAF", "pop_GBR_MAF"))
info(header(vcfHeaderExist)) <- rbind(info(header(vcfHeaderExist)), newInfoHeader)

vcfDataExist <- vcf
newInfoData <- DataFrame(
    MAF = seq_along(vcfHeaderExist),
    pop_GBR_MAF = rev(seq_along(vcfHeaderExist))
)
info(vcfDataExist) <- cbind(info(vcfDataExist), newInfoData)

# Signatures ----

test_that("addFrequencies supports all signatures",{

    # \alias{addFrequencies,ExpandedVCF,list,TVTBparam-method}
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

    # \alias{addFrequencies,ExpandedVCF,character,TVTBparam-method}
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

    # \alias{addFrequencies,ExpandedVCF,missing,TVTBparam-method}
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

# Overwrite INFO header fields ----

test_that(".checkFrequencyInfo overwrites INFO header fields", {

    expect_error(
        addFrequencies(
            vcf = vcfHeaderExist, param = tparam)
    )

    # expect_error(
    #     addFrequencies(
    #         vcf = vcfHeaderExist,
    #         phenos = list(pop = "GBR"),
    #         param = tparam)
    # )

    expect_message(
        addFrequencies(
            vcf = vcfHeaderExist, param = tparam, force = TRUE)
    )

    expect_message(
        addFrequencies(
            vcf = vcfHeaderExist,
            phenos = list(pop = "GBR"),
            param = tparam, force = TRUE)
    )

})

# Overwrite INFO data fields ----

test_that(".checkFrequencyInfo overwrites INFO data fields", {

    expect_error(
        addFrequencies(
            vcf = vcfDataExist, param = tparam)
    )

    # expect_error(
    #     addFrequencies(
    #         vcf = vcfDataExist,
    #         phenos = list(pop = "GBR"),
    #         param = tparam)
    # )

    # All the below produce a warning due to validity check of ExpandedVCF
    expect_warning(
        addFrequencies(
            vcf = vcfDataExist, param = tparam, force = TRUE)
    )

    expect_warning(
        addFrequencies(
            vcf = vcfDataExist,
            phenos = list(pop = "GBR"),
            param = tparam, force = TRUE)
    )

})
