context("addCountGenos")

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
vcf <- VariantAnnotation::expand(vcf)

# Arguments ----

test_that("addCountGenos supports all signatures",{

    ## Implicitely tested by higher-level methods
    # samples = "missing"
    expect_s4_class(
        addCountGenos(
            vcf = vcf, genos = het(tparam),
            key = "NHET",
            description = "Number of heterozygous genotypes"),
        "ExpandedVCF"
    )

    # samples = "numeric"
    expect_s4_class(
        addCountGenos(
            vcf = vcf, genos = c("0|1", "1|0"),
            key = "NHET",
            description = "Number of heterozygous genotypes",
            samples = 1:ncol(vcf)),
        "ExpandedVCF"
    )

    # samples = "character"
    expect_s4_class(
        addCountGenos(
            vcf = vcf, genos = c("0|1", "1|0"),
            key = "NHET",
            description = "Number of heterozygous genotypes",
            samples = colnames(geno(vcf)[["GT"]])),
        "ExpandedVCF"
    )

})

vcf_NHET <- addCountGenos(
    vcf = vcf, genos = c("0|1", "1|0"),
    key = "NHET",
    description = "Number of heterozygous genotypes")

# Argument: force ----

test_that("force=FALSE throws an error if the field exists",{

    expect_error(
        addCountGenos(
            vcf = vcf_NHET, genos = c("0|1", "1|0"),
            key = "NHET",
            description = "Number of heterozygous genotypes")
    )

})

test_that("force=TRUE messages that the field will be updated",{

    expect_message(
        addCountGenos(
            vcf = vcf_NHET, genos = c("0|1", "1|0"),
            key = "NHET",
            description = "Number of heterozygous genotypes",
            samples = colnames(geno(vcf)[["GT"]]),
            force = TRUE)
    )

})
