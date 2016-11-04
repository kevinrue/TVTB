context("addPhenoLevelFrequencies")

# Settings ----

# Genomic region
bedRegions <- GRanges(
    seqnames = "15",
    ranges = IRanges(start = 48420E3, end = 48421E3))

# VCF file
extdata <- file.path(system.file(package = "tSVE"), "extdata")
vcfFile <- file.path(extdata, "chr15.phase3_integrated.vcf.gz")
tabixVcf <- TabixFile(file = vcfFile)

# Good and bad phenotype files
phenoFile <- file.path(extdata, "integrated_samples.txt")
phenotypes <- DataFrame(read.table(
    file = phenoFile, header = TRUE, row.names = 1))
# Subset phenotypes to test with a small number of samples
samplePhenotypes <- subset(phenotypes, pop == "GBR")

# Import variants
svp <- ScanVcfParam(
    fixed = "ALT",
    info = "CSQ",
    geno = "GT",
    samples = rownames(samplePhenotypes),
    which = bedRegions)
vcf <- readVcf(file = tabixVcf, param = svp)
# Separate multi-allelic records into bi-allelic records
eVcf <- expand(x = vcf, row.names = TRUE)
# Disambiguate row.names from multi-allelic records
rownames(eVcf) <- paste(rownames(eVcf), mcols(eVcf)[,"ALT"], sep = "_")
# Add some phenotypes information necessary for the demo
colData(eVcf) <- samplePhenotypes

# Define genotypes ----
tparam <- tSVEParam(genos = list(
    REF = c("0|0"),
    HET = c("0|1","1|0"),
    ALT = c("1|1")))

# Signatures ----

test_that("addFrequencies supports all signatures",{

    # \alias{addPhenoLevelFrequencies,ExpandedVCF,tSVEParam-method}
    expect_s4_class(
        addPhenoLevelFrequencies(
            vcf = eVcf, pheno = "pop", level = "GBR", param = tparam),
        "ExpandedVCF"
    )

    # \alias{addPhenoLevelFrequencies,ExpandedVCF,missing-method}
    expect_s4_class(
        addPhenoLevelFrequencies(
            vcf = eVcf, pheno = "pop", level = "GBR",
            ref = "0|0", het = c("0|1", "1|0"), alt = "1|1"),
        "ExpandedVCF"
    )

})
