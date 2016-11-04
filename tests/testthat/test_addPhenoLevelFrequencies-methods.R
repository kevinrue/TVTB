context("addPhenoLevelFrequencies")

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

test_that("addPhenoLevelFrequencies supports all signatures",{

    # \alias{addPhenoLevelFrequencies,ExpandedVCF,TVTBparam-method}
    expect_s4_class(
        addPhenoLevelFrequencies(
            vcf = vcf, pheno = "pop", level = "GBR", param = tparam),
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
            vcf = vcfHeaderExist, pheno = "pop", level = "GBR", param = tparam)
    )

    expect_error(
        addPhenoLevelFrequencies(
            vcf = vcfDataExist, pheno = "pop", level = "GBR", param = tparam)
    )

})
