context("Use ensembldb package")

if(!requireNamespace("EnsDb.Hsapiens.v75")){
    stop("Package 'EnsDb.Hsapiens.v75' required for testing")
}

test_that("getEdb() returns appropriate values",{

    # Check that the function returns a value with proper input
    expect_s4_class(
        getEdb(x = "EnsDb.Hsapiens.v75"),
        "EnsDb")

})

test_that("TheFilter() returns appropriate values",{

    # Check that
    expect_s4_class(
        EnsDbFilter(type = "Genename", condition = "=", value = "SLC24A5"),
        "GenenameFilter"
    )

    # Check that
    expect_s4_class(
        EnsDbFilter(
            type = "Genename",
            condition = "=",
            value = "SLC24A5,ACTN1"),
        "GenenameFilter"
    )

    # Check that
    expect_s4_class(
        EnsDbFilter(
            type = "Genename",
            condition = "=",
            value = "SLC24A5 ACTN1"),
        "GenenameFilter"
    )

    # Check that
    expect_warning(
        EnsDbFilter(type = "Invalid", condition = "=", value = "SLC24A5")
    )

})


