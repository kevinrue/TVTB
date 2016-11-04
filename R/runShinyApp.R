# running the shiny web app.

tSVE <- function(...){
    if(requireNamespace("shiny", quietly=TRUE)){
        message("Starting the Shiny web app.")
        shiny::runApp(appDir = system.file("shinyApp", package = "TVTB"), ...)
    }else{
        stop("Package shiny not installed!")
    }
}
