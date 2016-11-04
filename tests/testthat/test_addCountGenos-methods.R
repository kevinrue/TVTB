context("addCountGenos")

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

# Arguments ----

test_that("addCountGenos supports all signatures",{

    ## Implicitely tested by higher-level methods
    # samples = "missing"
    expect_s4_class(
        addCountGenos(
            vcf = eVcf, genos = c("0|1", "1|0"),
            key = "NHET",
            description = "Number of heterozygous genotypes"),
        "ExpandedVCF"
    )

    # samples = "numeric"
    expect_s4_class(
        addCountGenos(
            vcf = eVcf, genos = c("0|1", "1|0"),
            key = "NHET",
            description = "Number of heterozygous genotypes",
            samples = 1:ncol(eVcf)),
        "ExpandedVCF"
    )

    # samples = "character"
    expect_s4_class(
        addCountGenos(
            vcf = eVcf, genos = c("0|1", "1|0"),
            key = "NHET",
            description = "Number of heterozygous genotypes",
            samples = colnames(geno(eVcf)[["GT"]])),
        "ExpandedVCF"
    )

})

eVcf_NHET <- addCountGenos(
    vcf = eVcf, genos = c("0|1", "1|0"),
    key = "NHET",
    description = "Number of heterozygous genotypes")

# Argument: force ----

test_that("force=FALSE throws an error if the field exists",{

    expect_error(
        addCountGenos(
            vcf = eVcf_NHET, genos = c("0|1", "1|0"),
            key = "NHET",
            description = "Number of heterozygous genotypes")
    )

})

test_that("force=TRUE messages that the field will be updated",{

    expect_message(
        addCountGenos(
            vcf = eVcf_NHET, genos = c("0|1", "1|0"),
            key = "NHET",
            description = "Number of heterozygous genotypes",
            samples = colnames(geno(eVcf)[["GT"]]),
            force = TRUE)
    )

})
