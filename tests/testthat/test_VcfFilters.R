context("filterExpandedVcf")

filterDataFrame <- DataFrame(
    minMaf = sample(x = c(TRUE, FALSE), size = 15, replace = TRUE),
    maxMaf = sample(x = c(TRUE, FALSE), size = 15, replace = TRUE),
    row.names = sample(x = LETTERS, size = 15)
)

print("filterDataFrame")
print(as.data.frame(filterDataFrame))

test_that(".filterExpandedVcf uses FilterRules correctly", {

    print("filterDataFrame")
    print(as.data.frame(filterDataFrame))

    filterRules <- .extractFilterRules(filterDataFrame = filterDataFrame)

    filterResults <- eval(expr = filterRules, envir = filterDataFrame)

    keepIdx <- which(filterResults)
})
