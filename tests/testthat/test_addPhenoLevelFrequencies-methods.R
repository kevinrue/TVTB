context("addPhenoLevelFrequencies")

# Settings ----

# VCF file
vcfFile <- system.file("extdata", "moderate.vcf", package = "TVTB")

# Phenotype file
phenoFile <- system.file("extdata", "moderate_pheno.txt", package = "TVTB")
phenotypes <- S4Vectors::DataFrame(read.table(phenoFile, TRUE, row.names = 1))

# TVTB parameters
tparam <- TVTBparam(Genotypes("0|0", c("0|1", "1|0"), "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(
    vcfFile, param = tparam, colData = phenotypes)
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)

# Create a VCF object with a pre-existing INFO key

vcfHeaderExist <- vcf
newInfoHeader <- S4Vectors::DataFrame(
    Number = rep(1, 2),
    Type = "Integer",
    Description = "Pre-existing INFO field",
    row.names = c("MAF", "pop_GBR_MAF"))
info(header(vcfHeaderExist)) <- rbind(
    info(header(vcfHeaderExist)), newInfoHeader)

vcfDataExist <- vcfHeaderExist
newInfoData <- S4Vectors::DataFrame(
    MAF = seq_along(vcfHeaderExist),
    pop_GBR_MAF = rev(seq_along(vcfHeaderExist))
)
info(vcfDataExist) <- cbind(info(vcfDataExist), newInfoData)


# Signatures ----

test_that("addPhenoLevelFrequencies supports all signatures",{

    # \alias{addPhenoLevelFrequencies,ExpandedVCF-method}
    expect_s4_class(
        addPhenoLevelFrequencies(vcf, "pop", "GBR"),
        "ExpandedVCF"
    )

})

# .checkInputsPLF ----

test_that(".checkInputsPLF catches invalid inputs", {

    # Invalid phenotype name
    expect_error(addPhenoLevelFrequencies(vcf, "missing", "GBR"))

    # Invalid phenotype level name
    expect_error(addPhenoLevelFrequencies(vcf, "pop", "missing"))

})

# Existing INFO header fields ----

test_that("Error thrown if INFO key cannot be overwritten", {

    expect_error(
        addPhenoLevelFrequencies(vcfHeaderExist, "pop", "GBR")
    )

    expect_error(
        addPhenoLevelFrequencies(vcfDataExist, "pop", "GBR")
    )

})

# Existing INFO data ----

test_that("Messages when overwriting INFO key", {

    # In both header and data
    expect_message(
        addPhenoLevelFrequencies(
            vcfDataExist, "pop", "GBR", force = TRUE
        )
    )

})
