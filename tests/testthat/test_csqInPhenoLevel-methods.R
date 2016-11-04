context("csqInPhenoLevel")

# Settings ----

# Genomic region
bedRegions <- GRanges(
    seqnames = "15",
    ranges = IRanges(start = 48413170, end = 48434757))

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

tparam <- tSVEParam(genos = list(c("0|0"), c("0|1","1|0"), c("1|1")))

# Signatures ----

test_that("csqInPhenoLevel() supports all signatures",{

    # ExpandedVCF, / implicitely tested by higher functions TODO:
    expect_is(
        csqInPhenoLevel(
            vcf = eVcf, phenoCol = "super_pop", level = "EUR",
            csqCol = "CADD_PHRED", param = tparam,
            het = c("0/1","1/0"),
            alt = "1/1",
            unique = FALSE, facet = NULL),
        "data.frame"
    )

    # ExpandedVCF, character
    expect_is(
        csqInPhenoLevel(
            vcf = eVcf, phenoCol = "super_pop", level = "EUR",
            csqCol = "CADD_PHRED",
            alts = unlist(carrier(tparam)),
            unique = FALSE, facet = NULL),
        "data.frame"
    )

})
