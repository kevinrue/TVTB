context("densityCsqByPhenotype")

# Settings ----

# Genomic region
regions <- GenomicRanges::GRanges(
    seqnames = "15",
    ranges = IRanges(start = 48413170, end = 48434757)) # all variants

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
    which = regions)
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

test_that("densityCsq* supports all signatures",{

    ## ByPhenotype()
    expect_is(
        densityCsqByPhenotype(
            vcf = eVcf,
            phenoCol = "gender",
            csqCol = "CADD_PHRED",
            param = tparam),
        "data.frame"
    )



    expect_is(
        densityCsqByPhenotype(
            vcf = eVcf,
            phenoCol = "gender",
            csqCol = "CADD_PHRED",
            alts = c("0|1", "1|0", "1|1")),
        "data.frame"
    )

    ## InPhenoLevel
    expect_is(
        densityCsqInPhenoLevel(
            level = "GBR",
            vcf = eVcf,
            phenoCol = "pop",
            csqCol = "CADD_PHRED",
            param = tparam),
        "data.frame"
    )

    ## Implicitely tested by *ByPhenotype
    expect_is(
        densityCsqInPhenoLevel(
            level = "GBR",
            vcf = eVcf,
            phenoCol = "pop",
            csqCol = "CADD_PHRED",
            param = tparam),
        "data.frame"
    )

})

# Arguments ----

test_that("plot & facet & popFreq argument work",{

    expect_s3_class(
        densityCsqByPhenotype(
            vcf = eVcf,
            phenoCol = "super_pop",
            csqCol = "AMR_MAF",
            param = tparam,
            facet = "Feature", plot = TRUE, popFreq = TRUE),
        c("gg", "ggplot")
    )

    expect_s3_class(
        densityCsqInPhenoLevel(
            level = "EUR",
            vcf = eVcf,
            phenoCol = "super_pop",
            csqCol = "AMR_MAF",
            param = tparam,
            facet = "Feature", plot = TRUE, popFreq = TRUE),
        c("gg", "ggplot")
    )

})
