# running the shiny web app.

tSVE <- function(...){ # nocov start
    if (requireNamespace("shiny", quietly=TRUE)){
        message("Starting the Shiny web app.")
        shiny::runApp(system.file("shinyApp", package = "TVTB"), ...)
    } else {
        stop("Package shiny not installed!")
    }
} # nocov end
