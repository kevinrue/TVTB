
scriptFolder <- file.path(system.file(package = "TVTB"), "shinyApp")

source(file.path(scriptFolder, "regions2vcf.R"))
source(file.path(scriptFolder, "fileParseFunctions.R"))
source(file.path(scriptFolder, "grangesParseFunctions.R"))
source(file.path(scriptFolder, "ensembldbFunctions.R"))
