
# Constructors ----

test_that("Constructor produce a valid object",{

    expect_s4_class(
        Genotypes("0/0", c("0/1", "1/0"), "1/1"),
        "Genotypes"
    )

})

genotypes <- Genotypes("0/0", c("0/1", "1/0"), "1/1")

# Validity checks ----

test_that("Genotypes must have lengths > 0", {

    expect_error(Genotypes(character(), c("0/1", "1/0"), "1/1"))

    expect_error(Genotypes("0/0", character(), "1/1"))

    expect_error(Genotypes("0/0", c("0/1", "1/0"), character()))

    expect_error(Genotypes("0/0", c("0/1", "1/0"), "1/1", c("REF", "HET")))

    expect_error(Genotypes("0/0", c("0/1", "1/0"), "1/1",c("REF", "HET", NA)))

})

# Accessors ----

test_that("Accessors return valid values", {

    expect_type(
        genos(genotypes),
        "character"
    )

    expect_type(
        ref(genotypes),
        "character"
    )

    expect_type(
        het(genotypes),
        "character"
    )

    expect_type(
        alt(genotypes),
        "character"
    )

    expect_type(
        carrier(genotypes),
        "character"
    )

    expect_type(
        suffix(genotypes),
        "character"
    )

    expect_identical(
        names(suffix(genotypes)),
        c("ref", "het", "alt")
    )

})

# Setters ----

test_that("Setters return valid values", {

    expect_type(
        ref(genotypes) <- "d",
        "character"
    )

    expect_type(
        het(genotypes) <- "e",
        "character"
    )

    expect_type(
        alt(genotypes) <- "f",
        "character"
    )

})

# Overlap ----

test_that("Overlapping genotypes are not allowed",{

    expect_error(Genotypes("0/0", c("0/1", "0/0"), "1/1"))

})

# show ----

test_that("show method return valid values",{

    expect_null(show(genotypes))

})
