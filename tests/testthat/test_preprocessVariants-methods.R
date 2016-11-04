context("preprocessVariants")

# Settings ----

# Genomic region
regions <- GRanges(
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
samplePhenotypes <- phenotypes[
    sample(x = 1:nrow(phenotypes), size = 100),]

# Bad phenotype file (not all phenoSamples are present in VCF)
phenoFilePartial <- file.path(system.file(
    package = "tSVE"), "badexamples", "pheno_PartialOverlap.txt")
phenotypesPartial <- DataFrame(read.table(
    file = phenoFilePartial, header = TRUE, row.names = 1))

# Bad phenotype file (no phenoSample is present in VCF)
phenoFileNoOverlap <- file.path(system.file(
    package = "tSVE"), "badexamples", "pheno_NoOverlap.txt")
phenotypesNoOverlap <- DataFrame(read.table(
    file = phenoFileNoOverlap, header = TRUE, row.names = 1))

# Actually, only the vep parameter is used at this step
tparam <- tSVEParam(
    ref = c("0|0"),
    het = c("0|1","1|0"),
    alt = c("1|1"),
    vep = "CSQ")

# Signatures ----

test_that("preprocessVariants supports all signatures",{

    ### file=TabixFile
    ## param=tSVEParam

    # file=TabixFile,param=tSVEParam,phenos=DataFrame
    expect_s4_class(
        preprocessVariants(
            file = tabixVcf, regions = regions,
            param = tparam, phenos = samplePhenotypes),
        "ExpandedVCF"
    )

    # file=TabixFile,param=tSVEParam,phenos=data.frame
    expect_s4_class(
        preprocessVariants(
            file = tabixVcf, regions = regions,
            param = tparam, phenos = as.data.frame(samplePhenotypes)),
        "ExpandedVCF"
    )

    # file=TabixFile,param=tSVEParam,phenos=character
    expect_s4_class(
        preprocessVariants(
            file = tabixVcf, regions = regions,
            param = tparam, phenos = phenoFile),
        "ExpandedVCF"
    )

    # file=TabixFile,param=tSVEParam,phenos=missing
    expect_s4_class(
        preprocessVariants(
            file = tabixVcf, regions = regions,
            param = tparam),
        "ExpandedVCF"
    )

    ## param=missing

    # file=TabixFile,param=missing,phenos=DataFrame
    expect_s4_class(
        preprocessVariants(
            file = tabixVcf, regions = regions,
            phenos = samplePhenotypes),
        "ExpandedVCF"
    )

    # file=TabixFile,param=missing,phenos=data.frame
    expect_s4_class(
        preprocessVariants(
            file = tabixVcf, regions = regions,
            phenos = as.data.frame(samplePhenotypes)),
        "ExpandedVCF"
    )

    # file=TabixFile,param=missing,phenos=character
    expect_s4_class(
        preprocessVariants(
            file = tabixVcf, regions = regions,
            phenos = phenoFile),
        "ExpandedVCF"
    )

    # file=TabixFile,param=missing,phenos=missing
    expect_s4_class(
        preprocessVariants(
            file = tabixVcf, regions = regions),
        "ExpandedVCF"
    )

    ### file=character
    ## param=tSVEParam

    # file=character,param=tSVEParam,phenos=DataFrame
    expect_s4_class(
        preprocessVariants(
            file = path(tabixVcf), regions = regions,
            param = tparam, phenos = samplePhenotypes),
        "ExpandedVCF"
    )

    # file=character,param=tSVEParam,phenos=data.frame
    expect_s4_class(
        preprocessVariants(
            file = path(tabixVcf), regions = regions,
            param = tparam, phenos = as.data.frame(samplePhenotypes)),
        "ExpandedVCF"
    )

    # file=character,param=tSVEParam,phenos=character
    expect_s4_class(
        preprocessVariants(
            file = path(tabixVcf), regions = regions,
            param = tparam, phenos = phenoFile),
        "ExpandedVCF"
    )

    # file=character,param=tSVEParam,phenos=missing
    expect_s4_class(
        preprocessVariants(
            file = path(tabixVcf), regions = regions,
            param = tparam),
        "ExpandedVCF"
    )

    ## param=missing

    # file=character,param=missing,phenos=DataFrame
    expect_s4_class(
        preprocessVariants(
            file = path(tabixVcf), regions = regions,
            phenos = samplePhenotypes),
        "ExpandedVCF"
    )

    # file=character,param=missing,phenos=data.frame
    expect_s4_class(
        preprocessVariants(
            file = path(tabixVcf), regions = regions,
            phenos = as.data.frame(samplePhenotypes)),
        "ExpandedVCF"
    )

    # file=character,param=missing,phenos=character
    expect_s4_class(
        preprocessVariants(
            file = path(tabixVcf), regions = regions,
            phenos = phenoFile),
        "ExpandedVCF"
    )

    # file=character,param=missing,phenos=missing
    expect_s4_class(
        preprocessVariants(
            file = path(tabixVcf), regions = regions),
        "ExpandedVCF"
    )

})

# test_that("vep are mandatory in the absence of tSVEParam",{
#
#     ## param=tSVEParam
#
#     # file=TabixFile,phenos=DataFrame,param=missing
#     expect_error(
#         preprocessVariants(
#             file = tabixVcf, regions = regions, phenos = samplePhenotypes,
#             vep = "CSQ")
#     )
#
# })

# Invalid phenos arguments ----

test_that("invalid phenos= arguments are detected",{

    expect_error(
        preprocessVariants(
            file = tabixVcf, regions = regions,
            param = tparam,
            phenos = phenotypesPartial)
    )

    expect_error(
        preprocessVariants(
            file = tabixVcf, regions = regions,
            param = tparam,
            phenos = phenotypesNoOverlap)
    )

})

# Invalid csqField arguments ----

test_that("invalid vep arguments are detected",{

    expect_error(
        preprocessVariants(
            file = tabixVcf, regions = regions, vep = "ERROR",
            param = tparam)
    )

})
