context("TVTBparam")

# Settings ----

genotypes <- Genotypes("0|0", c("0|1", "1|0"), "1|1")

gr <- GenomicRanges::GRanges(
    seqnames = "15", ranges = IRanges::IRanges(
        start = 48413170,
        end = 48434757,
        names = "SLC24A5"))

grl <- GenomicRanges::GRangesList(gr)

# Constructors ----

test_that("Constructors produce a valid object",{

    expect_s4_class(
        TVTBparam(genotypes),
        "TVTBparam"
    )

    expect_s4_class(
        TVTBparam(),
        "TVTBparam"
    )

})

# Conflict of suffixes ----

test_that("suffixes cannot overlap",{

    expect_error(
        TVTBparam(genotypes, aaf = "REF")
    )

})

# Conflict vep/INFO ----

test_that("vep key must be in ScanVcfParam:INFO",{

    expect_error(
        TVTBparam(genotypes, vep = "CSQ", svp = ScanVcfParam(info = "ANN"))
    )

})

# Accessors ----

tparam <- TVTBparam(genotypes)

test_that("Accessors return valid values",{

    expect_s4_class(
        genos(tparam),
        "Genotypes"
    )

    expect_type(
        ref(tparam),
        "character"
    )

    expect_type(
        het(tparam),
        "character"
    )

    expect_type(
        alt(tparam),
        "character"
    )

    expect_type(
        carrier(tparam),
        "character"
    )

    expect_s4_class(
        ranges(tparam),
        "GRangesList"
    )

    expect_type(
        aaf(tparam),
        "character"
    )

    expect_type(
        maf(tparam),
        "character"
    )

    expect_type(
        vep(tparam),
        "character"
    )

    expect_s4_class(
        svp(tparam),
        "ScanVcfParam"
    )

    expect_type(
        suffix(tparam),
        "character"
    )

    expect_s4_class(
        bp(tparam),
        "BiocParallelParam"
    )

})

# Setters ----

test_that("Setters return valid values",{

    expect_s4_class(
        genos(tparam) <- Genotypes("a", "b", "c"),
        "Genotypes"
    )

    expect_type(
        ref(tparam) <- "0/0",
        "character"
    )

    expect_type(
        het(tparam) <- c("0/1", "1/0"),
        "character"
    )

    expect_type(
        alt(tparam) <- c("1/1"),
        "character"
    )

    expect_s4_class(
        ranges(tparam) <- grl,
        "GRangesList"
    )

    expect_type(
        aaf(tparam) <- "aaf",
        "character"
    )

    expect_type(
        maf(tparam) <- "maf",
        "character"
    )

    expect_type(
        maf(tparam) <- "VEP",
        "character"
    )

    expect_type(
        vep(tparam) <- "CSQ",
        "character"
    )

    expect_s4_class(
        svp(tparam) <- ScanVcfParam(fixed = c("REF", "ALT")),
        "ScanVcfParam"
    )

    expect_s4_class(
        bp(tparam) <- BiocParallel::SerialParam(),
        "BiocParallelParam"
    )

})

# Invalid setters inputs ----

test_that("Setters catch invalid inputs",{

    expect_warning(
        alt(tparam) <- NA_character_
    )

    expect_error(
        aaf(tparam) <- c("aaf", "maf")
    )

    expect_error(
        maf(tparam) <- c("aaf", "maf")
    )

    expect_error(
        vep(tparam) <- c("aaf", "maf", "vep")
    )

})

# show ----

test_that("show method return valid values",{

    expect_null(show(tparam))

})

# Coerce ----

test_that("TVTBparam can be coerced to ScanVcfParam",{

    expect_s4_class(
        as(tparam, "ScanVcfParam"),
        "ScanVcfParam"
    )

})
