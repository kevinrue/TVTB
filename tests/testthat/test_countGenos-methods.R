context("countGenos")

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
vcf <- VariantAnnotation::expand(vcf)

# countGenos() ---

test_that("countGenos() returns appropriate values",{

    expect_type(
        countGenos(
            x = vcf,
            genos = unlist(het(tparam), use.names = FALSE),
            pheno = "pop",
            level = "GBR"),
        "integer"
    )

    expect_type(
        countGenos(
            x = geno(vcf)[["GT"]],
            genos = unlist(het(tparam), use.names = FALSE)),
        "integer"
    )

})

# .checkPhenoLevel ----

test_that(".checkPhenoLevel() catches invalid inputs",{

    expect_error(
        countGenos(
            x = vcf,
            genos = unlist(het(tparam), use.names = FALSE),
            pheno = "missing",
            level = "GBR")
    )

    expect_error(
        countGenos(
            x = vcf,
            genos = unlist(het(tparam), use.names = FALSE),
            pheno = "pop",
            level = "missing")
    )

    expect_error(
        countGenos(
            x = vcf,
            genos = unlist(het(tparam), use.names = FALSE),
            level = "missing")
    )

})
