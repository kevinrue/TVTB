context("readVcf")

# Settings ----

# VCF file
vcfFile <- system.file("extdata", "moderate.vcf.gz", package = "TVTB")

# Phenotype file
phenoFile <- system.file("extdata", "moderate_pheno.txt", package = "TVTB")
phenotypes <- S4Vectors::DataFrame(read.table(phenoFile, TRUE, row.names = 1))

# TVTB parameters
tparam <- TVTBparam(Genotypes("0|0", c("0|1", "1|0"), "1|1"))

invalidVepParam <- TVTBparam(
    Genotypes("0|0", c("0|1", "1|0"), "1|1"),
    vep = "MULTI_ALLELIC")

# readVcf ----

test_that("method supports all signature", {

    # file=character
    expect_s4_class(
        readVcf(
            vcfFile, param = tparam, colData = phenotypes,
            autodetectGT = TRUE),
        "CollapsedVCF"
    )

    # file=TabixFile
    expect_s4_class(
        readVcf(
            Rsamtools::TabixFile(vcfFile), param = tparam,
            colData = phenotypes),
        "CollapsedVCF"
    )

})

test_that("vep field absent from VCF header throws an error", {
    vcfFile2 <- system.file(
        "extdata", "chr15.phase3_integrated.vcf.gz", package = "TVTB")
    expect_error(readVcf(vcfFile, param = invalidVepParam))

})

test_that("different samples in colData and ScanVcfParam throw an error", {
    VariantAnnotation::vcfSamples(svp(tparam)) <- LETTERS[1:10]
    expect_error(readVcf(vcfFile, param = tparam, colData = phenotypes))

})
