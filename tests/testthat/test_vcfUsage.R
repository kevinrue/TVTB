context("Pre-process variants from VCF")

# Settings ----

# Genomic region
bedRegions <- GenomicRanges::GRanges(
    seqnames = "15",
    ranges = IRanges::IRanges(start = 48420E3, end = 48421E3))

# VCF file
vcfFolder <- file.path(system.file(package = "tSVE"), "extdata")
vcfPattern <- "^chr%s\\..*\\.vcf\\.gz$"
vcfFilePattern <- gsub("%s", "15", vcfPattern)
vcfFile <- list.files(vcfFolder, vcfFilePattern, full.names = TRUE)
tabixVcf <- Rsamtools::TabixFile(file = vcfFile)

# Good phenotype file
# phenoFile <- file.path(
#     system.file(package = "tSVE"),
#     "extdata",
#     "integrated_samples.txt")
# phenotypes <- DataFrame(read.table(
#     file = phenoFile, header = TRUE, row.names = 1))

# chr2file() ----

test_that("chr2file() returns appropriate values",{

    # Check that the function returns a value with proper input
    expect_type(
        chr2file(
            chr = "15",
            pattern = vcfPattern,
            folder = vcfFolder),
        "character")

    # Even numeric chromosome should be given as character
    expect_error(chr2file(chr = 2, pattern = vcfPattern, folder = vcfFolder))

    # Pattern must contain %s
    expect_error(chr2file(
        chr = "15",
        pattern = "^chr1_.*\\.vcf\\.gz$",
        folder = vcfFolder))

    # Folder must exist
    expect_error(chr2file(
        chr = "15", pattern = vcfPattern, folder = "MissingFolder"))

    # No file matched in an existing folder
    expect_warning(chr2file(
        chr = "none", pattern = vcfPattern, folder = vcfFolder))
    expect_warning(chr2file(
        chr = "15",
        pattern = vcfPattern,
        folder = file.path(system.file(package = "tSVE"), "R")))

    # Multiple
    expect_error(chr2file(
        chr = "3",
        pattern = vcfPattern,
        folder = file.path(system.file(package = "tSVE"), "badexamples")))

})

# Import tests variants ----

# svp <- ScanVcfParam(
#     fixed = "ALT",
#     info = "CSQ",
#     geno = "GT",
#     samples = character(),
#     which = bedRegions)
# vcf <- readVcf(file = tabixVcf, param = svp)
# vcf <- expand(x = vcf, row.names = TRUE)
# # Separate multi-allelic records into bi-allelic records
# eVcf <- expand(x = vcf, row.names = TRUE)
# # Disambiguate row.names from multi-allelic records
# rownames(eVcf) <- paste(rownames(eVcf), mcols(eVcf)[,"ALT"], sep = "_")

# addOverallFrequencies() ---

# NOTE:
#     Useful tests to spot a problem, but functionality implicitely tested when
#     testing preprocessVariants(). No additional coverage

# test_that("addOverallFrequencies() catches invalid inputs",{
#
#     expect_s4_class(
#         addOverallFrequencies(
#             vcf = vcf, ref = c("0|0"), het = c("0|1", "1|0"),
#             alt = c("1|1")),
#         "ExpandedVCF"
#     )
#
#     vcfFreq <- addOverallFrequencies(
#         vcf = vcf, ref = c("0|0"), het = c("0|1", "1|0"),
#         alt = c("1|1"))
#
#     expect_type(
#         info(vcfFreq)[,"REF"],
#         "numeric"
#     )
#     expect_type(
#         info(vcfFreq)[,"HET"],
#         "numeric"
#     )
#     expect_type(
#         info(vcfFreq)[,"ALT"],
#         "numeric"
#     )
#     expect_type(
#         info(vcfFreq)[,"AAF"],
#         "numeric"
#     )
#     expect_type(
#         info(vcfFreq)[,"MAF"],
#         "numeric"
#     )
#
# })


# Add phenotype information ----

# Add some phenotypes information necessary for the demo
# colData(eVcf) <- phenotypes

# addPhenotypeFrequencies() ----

# test_that("addPhenotypeFrequencies() returns appropriate values",{
#
#     expect_s4_class(
#         addPhenotypeFrequencies(
#             vcf = eVcf, field = "super_pop",
#             ref = c("0|0"), het = c("0|1"), alt = c("1|1")),
#         "ExpandedVCF"
#     )
#
# })
