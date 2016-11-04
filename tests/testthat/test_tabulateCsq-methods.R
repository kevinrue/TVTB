context("densityCsqByPhenotype")

# Settings ----

# Genomic region
bedRegions <- GenomicRanges::GRanges(
    seqnames = "15",
    ranges = IRanges::IRanges(start = 48413170, end = 48434757)) # all variants

# VCF file
extdata <- file.path(system.file(package = "tSVE"), "extdata")
vcfFile <- file.path(extdata, "chr15.phase3_integrated.vcf.gz")
tabixVcf <- Rsamtools::TabixFile(file = vcfFile)

# Good and bad phenotype files
phenoFile <- file.path(extdata, "integrated_samples.txt")
phenotypes <- S4Vectors::DataFrame(read.table(
    file = phenoFile, header = TRUE, row.names = 1))
# Subset phenotypes to test with a small number of samples
samplePhenotypes <- subset(phenotypes, pop == "GBR")
samplePhenotypes <- droplevels(samplePhenotypes)

# Import variants
svp <- VariantAnnotation::ScanVcfParam(
    fixed = "ALT",
    info = "CSQ",
    geno = "GT",
    samples = rownames(samplePhenotypes),
    which = bedRegions)
vcf <- VariantAnnotation::readVcf(file = tabixVcf, param = svp)
# Separate multi-allelic records into bi-allelic records
eVcf <- VariantAnnotation::expand(x = vcf, row.names = TRUE)
# Disambiguate row.names from multi-allelic records
rownames(eVcf) <- paste(rownames(eVcf), mcols(eVcf)[,"ALT"], sep = "_")
# Add some phenotypes information necessary for the demo
SummarizedExperiment::colData(eVcf) <- samplePhenotypes

# Define genotypes ----
tparam <- tSVEParam(genos = list(
    REF = c("0|0"),
    HET = c("0|1","1|0"),
    ALT = c("1|1")))

# Signatures ----

test_that("tabulateCsq* supports all signatures",{

    # ExpandedVCF, tSVEParam
    expect_is(
        tabulateCsqByPhenotype(
            vcf = eVcf,
            phenoCol = "super_pop",
            csqCol = "Consequence",
            param = tparam,
            unique = TRUE),
        "data.frame"
    )

    expect_is(
        tabulateCsqByPhenotype(
            vcf = eVcf,
            phenoCol = "super_pop",
            csqCol = "Consequence",
            alts = c("0|1", "1|0", "1|1"),
            unique = TRUE),
        "data.frame"
    )

    ## Implicitely tested by *ByPhenotype
    expect_is(
        tabulateCsqInPhenoLevel(
            level = "EUR",
            vcf = eVcf,
            phenoCol = "super_pop",
            csqCol = "Consequence",
            param = tparam),
        "data.frame"
    )

    expect_is(
        tabulateCsqInPhenoLevel(
            level = "EUR",
            vcf = eVcf,
            phenoCol = "super_pop",
            csqCol = "Consequence",
            alts = c("0|1", "1|0", "1|1")),
        "data.frame"
    )

})

# Arguments ----

test_that("plot & facet & popFreq argument work",{

    expect_s3_class(
        tabulateCsqByPhenotype(
            vcf = eVcf, phenoCol = "super_pop",
            csqCol = "Consequence", param = tparam,
            facet = "Feature", plot = TRUE, percentage = TRUE) +
            theme(legend.text = element_text(size = rel(.5))),
        c("gg", "ggplot")
    )

    expect_s3_class(
        tabulateCsqInPhenoLevel(
            level = "EUR", vcf = eVcf,
            phenoCol = "super_pop",
            csqCol = "Consequence", param = tparam,
            facet = "Feature", plot = TRUE, percentage = TRUE) +
            theme(legend.text = element_text(size = rel(.5))),
        c("gg", "ggplot")
    )

    ## For 1% of extra coverage: plot=TRUE, percentage=FALSE,
    expect_s3_class(
        tabulateCsqByPhenotype(
            vcf = eVcf, phenoCol = "super_pop",
            csqCol = "Consequence", param = tparam,
            unique = FALSE, facet = "Feature",
            plot = TRUE) +
            theme(legend.text = element_text(size = rel(.5))),
        c("gg", "ggplot")
    )

    expect_s3_class(
        tabulateCsqInPhenoLevel(
            level = "EUR", vcf = eVcf,
            phenoCol = "super_pop",
            csqCol = "Consequence", param = tparam,
            unique = FALSE, facet = "Feature", plot = TRUE) +
            theme(legend.text = element_text(size = rel(.5))),
        c("gg", "ggplot")
    )
})

