context("countGenos")

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
