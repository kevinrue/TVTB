context("preprocessVariants")

# Settings ----

# Genomic region
gr <- GenomicRanges::GRanges(
    seqnames = "15",
    ranges = IRanges::IRanges(start = 48420E3, end = 48421E3))

# tSVEParam
tparam <- tSVEParam(
    ref = c("0|0"),
    het = c("0|1","1|0"),
    alt = c("1|1"))
# ... with ranges
tparam.gr <- tSVEParam
ranges(tparam) <- gr

# VCF file
extdata <- file.path(system.file(package = "tSVE"), "extdata")
vcfFile <- file.path(extdata, "moderate.vcf")
vcfgzFile <- paste(vcfFile, "gz", sep = ".")
tabixVcf <- Rsamtools::TabixFile(file = vcfgzFile)

# Good and bad phenotype files
phenoFile <- file.path(extdata, "moderate_pheno.txt")
phenotypes <- S4Vectors::DataFrame(read.table(
    file = phenoFile, header = TRUE, row.names = 1))

# Bad phenotype file (not all phenoSamples are present in VCF)
phenoFilePartial <- file.path(system.file(
    package = "tSVE"), "badexamples", "pheno_PartialOverlap.txt")
phenotypesPartial <- S4Vectors::DataFrame(read.table(
    file = phenoFilePartial, header = TRUE, row.names = 1))

# Bad phenotype file (no phenoSample is present in VCF)
phenoFileNoOverlap <- file.path(system.file(
    package = "tSVE"), "badexamples", "pheno_NoOverlap.txt")
phenotypesNoOverlap <- S4Vectors::DataFrame(read.table(
    file = phenoFileNoOverlap, header = TRUE, row.names = 1))

# Signatures ----

test_that("preprocessVariants supports all signatures",{

    ### file=TabixFile
    ## param=tSVEParam

    # file=TabixFile,param=tSVEParam,phenos=DataFrame
    expect_s4_class(
        preprocessVariants(
            file = tabixVcf,
            param = tparam, phenos = phenotypes),
        "ExpandedVCF"
    )

    # file=TabixFile,param=tSVEParam,phenos=data.frame
    expect_s4_class(
        preprocessVariants(
            file = tabixVcf,
            param = tparam, phenos = as.data.frame(phenotypes)),
        "ExpandedVCF"
    )

    # file=TabixFile,param=tSVEParam,phenos=character
    expect_s4_class(
        preprocessVariants(
            file = tabixVcf,
            param = tparam, phenos = phenoFile),
        "ExpandedVCF"
    )

    # file=TabixFile,param=tSVEParam,phenos=missing
    expect_s4_class(
        preprocessVariants(
            file = tabixVcf,
            param = tparam),
        "ExpandedVCF"
    )

    ## param=missing

    # file=TabixFile,param=missing,phenos=DataFrame
    expect_s4_class(
        preprocessVariants(
            file = tabixVcf,
            phenos = phenotypes),
        "ExpandedVCF"
    )

    # file=TabixFile,param=missing,phenos=data.frame
    expect_s4_class(
        preprocessVariants(
            file = tabixVcf,
            phenos = as.data.frame(phenotypes)),
        "ExpandedVCF"
    )

    # file=TabixFile,param=missing,phenos=character
    expect_s4_class(
        preprocessVariants(
            file = tabixVcf,
            phenos = phenoFile),
        "ExpandedVCF"
    )

    # file=TabixFile,param=missing,phenos=missing
    expect_s4_class(
        preprocessVariants(
            file = tabixVcf),
        "ExpandedVCF"
    )

    ### file=character
    ## param=tSVEParam

    # file=character,param=tSVEParam,phenos=DataFrame
    expect_s4_class(
        preprocessVariants(
            file = path(tabixVcf),
            param = tparam, phenos = phenotypes),
        "ExpandedVCF"
    )

    # file=character,param=tSVEParam,phenos=data.frame
    expect_s4_class(
        preprocessVariants(
            file = path(tabixVcf),
            param = tparam, phenos = as.data.frame(phenotypes)),
        "ExpandedVCF"
    )

    # file=character,param=tSVEParam,phenos=character
    expect_s4_class(
        preprocessVariants(
            file = path(tabixVcf),
            param = tparam, phenos = phenoFile),
        "ExpandedVCF"
    )

    # file=character,param=tSVEParam,phenos=missing
    expect_s4_class(
        preprocessVariants(
            file = path(tabixVcf),
            param = tparam),
        "ExpandedVCF"
    )

    ## param=missing

    # file=character,param=missing,phenos=DataFrame
    expect_s4_class(
        preprocessVariants(
            file = path(tabixVcf),
            phenos = phenotypes),
        "ExpandedVCF"
    )

    # file=character,param=missing,phenos=data.frame
    expect_s4_class(
        preprocessVariants(
            file = path(tabixVcf),
            phenos = as.data.frame(phenotypes)),
        "ExpandedVCF"
    )

    # file=character,param=missing,phenos=character
    expect_s4_class(
        preprocessVariants(
            file = path(tabixVcf),
            phenos = phenoFile),
        "ExpandedVCF"
    )

    # file=character,param=missing,phenos=missing
    expect_s4_class(
        preprocessVariants(
            file = path(tabixVcf)),
        "ExpandedVCF"
    )

})

# Invalid phenos arguments ----

test_that("invalid phenos= arguments are detected",{

    expect_error(
        preprocessVariants(
            file = tabixVcf,
            param = tparam,
            phenos = phenotypesPartial)
    )

    expect_error(
        preprocessVariants(
            file = tabixVcf,
            param = tparam,
            phenos = phenotypesNoOverlap)
    )

})

# Invalid csqField arguments ----

test_that("invalid vep arguments are detected",{

    expect_error(
        preprocessVariants(
            file = tabixVcf, vep = "ERROR",
            param = tparam)
    )

})

# file=character;ranges=GRanges(0) is incompatible ----

test_that("ranges cannot be used if VCF is not indexed",{

    expect_error(
        preprocessVariants(
            file = path(tabixVcf),
            param = tparam.gr, phenos = phenotypes)
    )

})
