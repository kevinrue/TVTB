context("addFrequencies")

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

# Overwrite INFO fields ----

test_that(".checkFrequencyInfo overwrites INFO fields", {

    expect_error(
        expect_s4_class(
            addFrequencies(
                vcf = vcfInfoExist, param = tparam),
            "ExpandedVCF"
        )
    )

    expect_message(
        expect_s4_class(
            addFrequencies(
                vcf = vcfInfoExist, param = tparam, force = TRUE),
            "ExpandedVCF"
        )
    )

    expect_error(
        expect_s4_class(
            addFrequencies(
                vcf = vcfInfoExist,
                phenos = list(pop = "GBR"),
                param = tparam),
            "ExpandedVCF"
        )
    )

    expect_message(
        expect_s4_class(
            addFrequencies(
                vcf = vcfInfoExist,
                phenos = list(pop = "GBR"),
                param = tparam, force = TRUE),
            "ExpandedVCF"
        )
    )



})
