context("addFrequencies")

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
vcf <- addFrequencies(vcf, "super_pop")

# Signatures ----

test_that("all signatures work to completion", {

    expect_type(
        plotInfo(
            vcf, "MAF",
            range(granges(vcf)),
            EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
            "super_pop"
        ),
        "list"
    )

})

# Helper .findInfoMetricsColumns ----

test_that("invalid metric/phenotype combination is detected", {

    expect_error(TVTB:::.findInfoMetricColumns(vcf, "INVALID_(.*)_invalid"))

})

# Helper: .geneRegionTrackFromEnsDb ----

test_that("invalid metric/phenotype combination is detected", {

    expect_s4_class(
        TVTB:::.geneRegionTrackFromEnsDb(
            EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75, range(granges(vcf))
        ),
        "GeneRegionTrack"
    )

})
