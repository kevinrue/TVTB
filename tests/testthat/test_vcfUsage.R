context("Pre-process variants from VCF")

# Settings ----

# VCF file
vcfFolder <- file.path(system.file(package = "TVTB"), "extdata")
vcfPattern <- "^chr%s\\..*\\.vcf\\.gz$"
vcfFilePattern <- gsub("%s", "15", vcfPattern)
vcfFile <- list.files(vcfFolder, vcfFilePattern, full.names = TRUE)
tabixVcf <- Rsamtools::TabixFile(file = vcfFile)

# Good phenotype file
# phenoFile <- file.path(
#     system.file(package = "TVTB"),
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
    expect_type(chr2file(
        chr = 15,
        pattern = vcfPattern,
        folder = vcfFolder),
        "character")

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
        folder = file.path(system.file(package = "TVTB"), "R")))

    # Multiple
    expect_warning(chr2file(
        chr = "3",
        pattern = vcfPattern,
        folder = file.path(system.file(package = "TVTB"), "badexamples")))

})
