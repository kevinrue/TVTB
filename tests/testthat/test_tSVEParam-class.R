context("tSVEParam")

# Settings ----

genotypes <- list(
    REF = c("0|0"),
    HET = c("0|1", "1|0"),
    ALT = c("1|1")
)

# Constructors ----

test_that("Constructors produce a valid object",{

    expect_s4_class(
        tSVEParam(genos = genotypes),
        "tSVEParam"
    )

    expect_s4_class(
        tSVEParam(
            ref = genotypes[["REF"]],
            het = genotypes[["HET"]],
            alt = genotypes[["ALT"]]),
        "tSVEParam"
    )

})

# Accessors ----

tparam <- new("tSVEParam", genos = genotypes)

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
        names(genos(tparam)) <- LETTERS[1:3],
        "character"
    )

    expect_type(
        hRef(tparam) <- list(ref = "0/0"),
        "list"
    )

    expect_type(
        het(tparam) <- list(het = c("0/1", "1/0")),
        "list"
    )

    expect_type(
        hAlt(tparam) <- list(alt = "1/1"),
        "list"
    )

    expect_type(
        carrier(tparam) <- list(
            het = c("0/1", "1/0"),
            alt = "1/1"),
        "list"
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

# Override ----

test_that("Override function return valid values",{

    expect_s4_class(
        .override.tSVEParam(
            param = tparam,
            ref = "0/0",
            het = "0/1",
            alt = "1/1",
            aaf = "Aaf",
            maf = "Maf",
            vep = "ANN",
            bp = MulticoreParam(workers = 2)),
        "tSVEParam"
    )

})
