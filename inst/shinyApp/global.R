
stopifnot(
  requireNamespace("BiocParallel"),
  requireNamespace("dplyr"),
  requireNamespace("DT"),
  requireNamespace("ensembldb"),
  requireNamespace("ensemblVEP"),
  requireNamespace("ggplot2"),
  requireNamespace("reshape2"),
  requireNamespace("Rsamtools"),
  requireNamespace("TVTB"),
  requireNamespace("VariantAnnotation"),
  requireNamespace("limma"),
  requireNamespace("rtracklayer")
)

# Display settings ----

# TODO: plot currently disabled
# .vepCountBarPlotHeight <- "500px"
# .vepDensityPlotHeight <- "500px"



# General settings --------------------------------------------------------

# Fetch the environment defined in the package
.tSVE <- TVTB:::.tSVE

# R session variables -----------------------------------------------------

.originalOptions <- options()
options("width" = get("options.width", .tSVE))

# Parallel settings -------------------------------------------------------

.PS <- list(
  default.bpClass = "SerialParam",
  choices.bpClass = names(BiocParallel::registered()),
  default.bpCores = 1,
  default.bpType = c(),
  choices.bpType = c()
)

# Messages ----------------------------------------------------------------

# Alphabetical order (ignoring ^default, ^choices, ^all, ...)
.Msgs <- list(
  # VCF
  importVariants = "Please import/refresh variants.",
  filteredVcfNULL = "Please wait while variants are being filtered...",
  singleVcf = "Please select a VCF file.",

  # GRanges
  noGenomicRanges = "No genomic range defined. All variants considered.",

  # Phenotypes
  phenoFile = "No phenotype file provided. All samples considered.",
  phenoInvalid = "Invalid phenotype data. See console.",
  colDataEmptyNoFile = "No phenotype information available.",
  colDataEmptyImported = paste(
    "Phenotypes imported, but not attached to variants.",
    "Please import/refresh variants to attach phenotypes."),

  # Specials
  fileChooseCancelled = "file choice cancelled.",
  vepKeyNotFound = "VEP key not found in INFO slot"
)

# Tracking messages -------------------------------------------------------

# Progress bars
.Tracking = list(
  # Calculations
  preprocessing = "Pre-processing data",
  calculate = "Crunching data",
  postprocessing = "Post-processing data",

  # Plotting
  ggplot = "Assembling plot (ggplot2)",
  render = "Rendering",

  # INFO frequencies
  rmFreqPhenoLevel = "Dropping phenotype-level frequencies",
  addFreqPhenoLevel = "Calculating phenotype-level frequencies",
  rmFreqOverall = "Dropping overall frequencies",
  addFreqOverall = "Calculating overall frequencies",

  # Imports
  singleVcf = "Importing from single VCF",
  multiVcfs = "Importing from multiple VCF files"
)

# EnsDb packages ----------------------------------------------------------

## list all packages...
.EnsDbPacks <- grep("^EnsDb", rownames(installed.packages()), value = TRUE)
