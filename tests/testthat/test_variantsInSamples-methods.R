context("variantsInSamples")

# Settings ----

# Genomic region
bedRegions <- GenomicRanges::GRanges(
    seqnames = "15",
    ranges = IRanges::IRanges(start = 48420E3, end = 48421E3))

# VCF file
extdata <- file.path(system.file(package = "tSVE"), "extdata")
vcfFile <- file.path(extdata, "chr15.phase3_integrated.vcf.gz")
tabixVcf <- Rsamtools::TabixFile(file = vcfFile)

# Good and bad phenotype files
phenoFile <- file.path(extdata, "integrated_samples.txt")
phenotypes <- S4Vectors::DataFrame(read.table(
    file = phenoFile, header = TRUE, row.names = 1))
# Subset phenotypes to test with a small number of samples
samplePhenotypes <- phenotypes[
    sample(x = 1:nrow(phenotypes), size = 100),]

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

sampleIdx <- 1:(ncol(eVcf)/10)

tparam <- tSVEParam(genos = list(c("0|0"), c("0|1","1|0"), c("1|1")))

# Signatures ----

test_that("variantsInSamples() supports all signatures",{

    ## ExpandedVCF,numeric,tSVEParam
    expect_type(
        variantsInSamples(
            vcf = eVcf,
            samples = sampleIdx,
            param = tparam,
            unique = FALSE),
        "integer"
    )

    # # ExpandedVCF,numeric,missing / tested by higher functions TODO: change
    expect_type(
        variantsInSamples(
            vcf = eVcf,
            samples = sampleIdx,
            alts = c("0|1", "1|0", "1|1"),
            unique = FALSE),
        "integer"
    )

    # ExpandedVCF,character,tSVEParam
    expect_type(
        variantsInSamples(
            vcf = eVcf,
            samples = colnames(eVcf)[sampleIdx],
            param = tparam,
            unique = FALSE),
        "integer"
    )

    # ExpandedVCF,character,missing
    expect_type(
        variantsInSamples(
            vcf = eVcf,
            samples = colnames(eVcf)[sampleIdx],
            alts = c("0|1", "1|0", "1|1"),
            unique = FALSE),
        "integer"
    )

})

# Argument: unique ----

test_that("unique=TRUE is supported",{

    # FALSE / tested above
    # expect_type(
    #     variantsInSamples(
    #         vcf = eVcf,
    #         samples = sampleIdx,
    #         alts = c("0|1", "1|0", "1|1"),
    #         unique = FALSE),
    #     "integer"
    # )

    # TRUE
    expect_type(
        variantsInSamples(
            vcf = eVcf,
            samples = sampleIdx,
            alts = c("0|1", "1|0", "1|1"),
            unique = TRUE),
        "integer"
    )

})

# Argument: unique ----

test_that("two or more alternate genotypes are required",{

    expect_error(
        variantsInSamples(
            vcf = eVcf,
            samples = sampleIdx,
            alts = c("1|1"),
            unique = FALSE)
    )

    expect_error(
        variantsInSamples(
            vcf = eVcf,
            samples = colnames(eVcf)[sampleIdx],
            alts = c("1|1"),
            unique = FALSE)
    )

})
