context("TVTBparam")

# Settings ----

genotypes <- list(
    REF = c("0|0"),
    HET = c("0|1", "1|0"),
    ALT = c("1|1")
)

genotypesNoName <- genotypes
names(genotypesNoName) <- NULL

genotypesPartiallyNamed <- genotypes
names(genotypesPartiallyNamed)[2] <- ""

genotypesOverlapping <- genotypes
genotypesOverlapping[[2]] <- genotypesOverlapping[[1]]

gr <- GenomicRanges::GRanges(
    seqnames = "15", ranges = IRanges(
        start = 48413170,
        end = 48434757,
        names = "SLC24A5"))

# Constructors ----

test_that("Constructors produce a valid object",{

    expect_s4_class(
        TVTBparam(genos = genotypes),
        "TVTBparam"
    )

    expect_s4_class(
        TVTBparam(
            ref = genotypes[["REF"]],
            het = genotypes[["HET"]],
            alt = genotypes[["ALT"]]),
        "TVTBparam"
    )

})

test_that("Constructors adds default genotype labels if missing",{

    expect_s4_class(
        TVTBparam(genos = genotypesNoName),
        "TVTBparam"
    )

})

# Invalid constructor inputs ----

test_that("Three genotypes are required",{

    expect_error(
        TVTBparam(genos = genotypes[1:2])
    )

})

test_that("Partially named genotypes are not allowed",{

    expect_error(
        TVTBparam(genos = genotypesPartiallyNamed)
    )

})

test_that("Overlapping genotypes are not allowed",{

    expect_error(
        TVTBparam(genos = genotypesOverlapping)
    )

})

# Accessors ----

tparam <- new("TVTBparam", genos = genotypes)

test_that("Accessors return valid values",{

    expect_type(
        genos(tparam),
        "list"
    )

    expect_type(
        sapply(genos(tparam), "class"),
        "character"
    )

    expect_type(
        hRef(tparam),
        "list"
    )

    expect_type(
        het(tparam),
        "list"
    )

    expect_type(
        hAlt(tparam),
        "list"
    )

    expect_type(
        carrier(tparam),
        "list"
    )

    expect_s4_class(
        ranges(tparam),
        "GRanges"
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

})

# Setters ----

test_that("Setters return valid values",{

    expect_type(
        genos(tparam) <- genotypes,
        "list"
    )

    expect_type(
        genos(tparam) <- genotypesNoName,
        "list"
    )

    expect_type(
        names(genos(tparam)) <- LETTERS[1:3],
        "character"
    )

    expect_type(
        hRef(tparam) <- list(ref = "0/0"),
        "list"
    )

    expect_type(
        hRef(tparam) <- "0/0",
        "character"
    )

    expect_type(
        het(tparam) <- list(het = c("0/1", "1/0")),
        "list"
    )

    expect_type(
        het(tparam) <- c("0/1", "1/0"),
        "character"
    )

    expect_type(
        hAlt(tparam) <- list(alt = "1/1"),
        "list"
    )

    expect_type(
        hAlt(tparam) <- c("1/1"),
        "character"
    )

    expect_type(
        carrier(tparam) <- list(
            het = c("0/1", "1/0"),
            alt = "1/1"),
        "list"
    )

    expect_s4_class(
        ranges(tparam) <- gr,
        "GRanges"
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

})

# Invalid setters inputs ----

test_that("Setters catch invalid inputs",{

    expect_error(
        genos(tparam) <- genotypesPartiallyNamed
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

# Override ----

test_that("Override method return valid values",{

    expect_s4_class(
        .override.TVTBparam(
            param = tparam,
            ref = "0/0",
            het = "0/1",
            alt = "1/1",
            ranges = gr,
            aaf = "Aaf",
            maf = "Maf",
            vep = "ANN",
            bp = MulticoreParam(workers = 2)),
        "TVTBparam"
    )

})

# show ----

test_that("show method return valid values",{

    expect_null(show(tparam))

})
