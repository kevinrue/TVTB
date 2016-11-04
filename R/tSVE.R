
# Create a locked environment to store the tSVE parameters
# Those parameters are set by the argument of the tSVE() method
.tSVE <- new.env(parent=emptyenv())
assign("refGT", "0|0", envir = .tSVE)
assign("hetGT", c("0|1", "1|2", "0|2", "1|0", "2|1", "2|0"), envir = .tSVE)
assign("altGT", c("1|1", "2|2"), envir = .tSVE)
assign("refSuffix", "REF", envir = .tSVE)
assign("hetSuffix", "HET", envir = .tSVE)
assign("altSuffix", "ALT", envir = .tSVE)
assign("mafSuffix", "MAf", envir = .tSVE)
assign("aafSuffix", "AAF", envir = .tSVE)
assign("vepKey", "CSQ", envir = .tSVE)
assign("genoHeatmap.height", "500px", envir = .tSVE)

lockEnvironment(.tSVE)

# nocov start
tSVE <- function(
    ...,
    refGT = "0|0",
    hetGT = c("0|1", "1|2", "0|2", "1|0", "2|1", "2|0"),
    altGT = c("1|1", "2|2"),
    vepKey = "CSQ",
    refSuffix = "REF", hetSuffix = "HET", altSuffix = "ALT",
    aafSuffix = "AAF", mafSuffix = "MAF",
    genoHeatmap.height = "500px"){
    if (requireNamespace("shiny", quietly=TRUE)){
        message("Setting environment ...")
        assign("refGT", refGT, envir = .tSVE)
        assign("hetGT", hetGT, envir = .tSVE)
        assign("altGT", altGT, envir = .tSVE)
        assign("refSuffix", refSuffix, envir = .tSVE)
        assign("hetSuffix", hetSuffix, envir = .tSVE)
        assign("altSuffix", altSuffix, envir = .tSVE)
        assign("aafSuffix", aafSuffix, envir = .tSVE)
        assign("mafSuffix", mafSuffix, envir = .tSVE)
        assign("vepKey", vepKey, envir = .tSVE)
        assign("genoHeatmap.height", genoHeatmap.height, envir = .tSVE)

        message("Starting the Shiny web app ...")
        shiny::runApp(system.file("shinyApp", package = "TVTB"), ...)
    } else {
        stop("Package shiny not installed!")
    }
}
# nocov end
