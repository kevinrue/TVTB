context("addCountGenos")

# Settings ----

# VCF file
vcfFile <- system.file("extdata", "moderate.vcf", package = "TVTB")

# TVTB parameters
tparam <- TVTBparam(Genotypes("0|0", c("0|1", "1|0"), "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(vcfFile, param = tparam)
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)

# Arguments ----

test_that("addCountGenos supports all signatures",{

    ## Implicitely tested by higher-level methods
    # samples = "missing"
    expect_s4_class(
        addCountGenos(
            vcf, het(tparam), "NHET", "Number of heterozygous genotypes"),
        "ExpandedVCF"
    )

    # samples = "numeric"
    expect_s4_class(
        addCountGenos(
            vcf, c("0|1", "1|0"),
            "NHET", "Number of heterozygous genotypes", samples = 1:ncol(vcf)),
        "ExpandedVCF"
    )

    # samples = "character"
    expect_s4_class(
        addCountGenos(
            vcf, c("0|1", "1|0"),
            "NHET", "Number of heterozygous genotypes",
            colnames(geno(vcf)[["GT"]])),
        "ExpandedVCF"
    )

})

vcf_NHET <- addCountGenos(
    vcf, c("0|1", "1|0"), "NHET", "Number of heterozygous genotypes")

# Argument: force ----

test_that("force=FALSE throws an error if the field exists",{

    expect_error(
        addCountGenos(
            vcf_NHET, c("0|1", "1|0"),
            "NHET", "Number of heterozygous genotypes")
    )

})

test_that("force=TRUE messages that the field will be updated",{

    expect_message(
        addCountGenos(
            vcf_NHET, c("0|1", "1|0"),
            "NHET", "Number of heterozygous genotypes",
            colnames(geno(vcf)[["GT"]]), TRUE)
    )

})
