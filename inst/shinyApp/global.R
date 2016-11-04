
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
    requireNamespace("VariantAnnotation")
)

message(
    "Please run ",
    "devtools::install_github(c(\"rstudio/shiny\", \"rstudio/DT\")) ",
    "to enable latest functionalities.")

# Seems like all the following don't need to be re-imported for the app
# library(GenomicRanges)
# library(reshape2)
# library(ensembldb)

# Display settings ----

vepCountBarplotHeight <- "500px"
vepDensityplotHeight <- "500px"

# R session variables -----------------------------------------------------

originalOptions <- options()
options("width" = 120)

# General settings --------------------------------------------------------

GS <- list(
    default.vep = "CSQ",
    default.refGenotypes = c("0|0"),
    default.hetGenotypes = c("0|1", "1|0", "0|2", "2|0", "1|2", "2|1"),
    default.altGenotypes = c("1|1", "2|2"),
    default.refSuffix = "REF", # Reference homozygote
    default.hetSuffix = "HET", # Heterozygote
    default.altSuffix = "ALT", # Alternate homozyote
    default.aafSuffix = "AAF", # Alternate allele frequency
    default.mafSuffix = "MAF", # Minor allele frequency
    choices.ensDbType = list(
        "Gene name" = "Genename"),
    default.ensDbType = "Genename",
    choices.ensDbFilters = c("=", "!=", "like", "in"),
    default.ensDbFilters = "=",
    all.genotypes = c(
        "0/0","0/1","1/0","1/1","./0","0/.","./1","1/.","./.",
        "0/2","2/0", "1/2", "2/1", "2/2", "./2", "2/.",
        "0|0","0|1","1|0", "1|1",".|0","0|.",".|1","1|.",
        "0|2","2|0", "1|2", "2|1", "2|2", ".|2", "2|.",
        "0", "1", "2","."),
    choices.phenoCols = c(),
    default.phenoCols = c(),
    choices.grangesInputMode = list(
        "BED file" = "bed",
        "UCSC browser" = "ucsc",
        "EnsDb package" = "EnsDb"
    ),
    default.grangesInputMode = c("bed"),
    choices.vcfCols = c(),
    default.vcfCols = c(),
    default.vcfFolder = system.file("extdata", package = "TVTB"),
    default.vcfPattern = "^chr%s\\..*\\.vcf\\.gz$",
    choices.vcfInputMode = list(
        "Single VCF" = "SingleVcf",
        "One per chromosome" = "OnePerChr"),
    default.vcfInputMode = "OnePerChr",
    vcfFilterClass.choices = list(
        "fixed" = "VcfFixedRules",
        "info" = "VcfInfoRules",
        "VEP" = "VcfVepRules"
    ),
    vcfFilterClass.default = "VcfFixedRules"
)

# Parallel settings -------------------------------------------------------

PS <- list(
    default.bpClass = "SerialParam",
    choices.bpClass = names(BiocParallel::registered()),
    default.bpCores = 1,
    default.bpType = c(),
    choices.bpType = c()
)

# Messages ----------------------------------------------------------------

# Alphabetical order (ignoring ^default, ^choices, ^all, ...)
Msgs <- list(
    # Mandatory inputs
    vepKey = "INFO field of VEP prediction must be defined.",
    refGenotypes = "Reference genotype(s) must be defined",
    hetGenotypes = "Heterozygote genotype(s) must be defined.",
    altGenotypes = "Alternative homozygote genotype(s) must be defined.",
    annotationPackage =
        "An EnsDb annotation package must be installed/selected.",
    invalidUcscRanges = "Invalid UCSC-type input.",
    vcfFolder = "VCF folder",
    vcfPattern = "VCF pattern",
    xAxisAngle = "X axis angle",
    xAxisHjust = "X axis horizontal justification",
    xAxisVjust = "X axis vertical justification",
    xAxisSize = "X axis font size",

    # INFO suffixes
    refSuffix = "Suffix for REF allele data",
    hetSuffix = "Suffix for HET allele data",
    altSuffix = "Suffix for ALT allele data",
    aafSuffix = "Suffix for ALT allele frequency",
    mafSuffix = "Suffix for minor allele frequency",

    # Tabulate VEP by phenotype (TVBP)
    vepTVBP = "VEP field tabulated",
    vepFacetKeyTVBP = "Faceting VEP field",
    stackedPercentageTVBP = "Stack percentage?",

    # Density VEP by phenotype (DVBP)
    vepDVBP = "VEP field tabulated",
    vepFacetKeyDVBP = "Faceting VEP field",

    # Optional inputs
    importVariants = "Please import/refresh variants.",
    noGenomicRanges = "No genomic range defined. All variants considered.",
    invalidGenomicRanges = "Invalid genomic range defined. See console.",
    phenotypes = "No phenotypes defined. All samples considered.",
    colDataEmptyOK = "No phenotype information available.",
    colDataEmptyImport = paste(
        "Phenotype information imported",
        "but not attached to VCF.",
        "Variants must be imported again to attach phenotype information."),
    singleVcf = "Please select a VCF file",

    # Calculated values (missing are likely bugs)
    filteredVcf = "Please filter variants.",
    filterVcfEmpty = "No variant left after filtering.",
    ensDbFilter = "Invalid EnsDb filter",
    genomeSeqinfo = "genomeSeqinfo bug?",
    genoSampleRanges.rows = "genoSampleRanges$rows bug?",
    genoSampleRanges.cols = "genoSampleRanges$cols bug?",
    genotypes = "genotypes bug?",
    legendText = "legendText bug?",
    mafRange = "mafRange bug?",
    phenotypes = "phenotypes bug?",
    # phenos = "phenos bug?",
    queryGenes = "queryGenes bug?",
    tfl = "tfl bug?",
    tparam = "tparam bug?",
    vcfContent = "vcfContent bug?",
    vcfFiles = "vcfFiles bug?",
    vcfHeader = "vcfHeader bug?",
    vcfPaths = "No VCF file found for requested chromosome. See console.",

    # Specials
    fileChooseCancelled = "file choice cancelled.",
    vepKeyNotFound = "VEP key not found in INFO slot"
)

# Tracking messages -------------------------------------------------------

# Progress bars
Tracking = list(
    # Calculations
    preprocessing = "Pre-processing data",
    calculate = "Crunching data...",
    postprocessing = "Post-processing data",

    # Plotting
    ggplot = "Assembling plot (ggplot2)",
    render = "Rendering.",

    # INFO frequencies
    rmFreqPhenoLevel = "Dropping phenotype-level frequencies",
    addFreqPhenoLevel = "Calculating phenotype-level frequencies",
    rmFreqOverall = "Dropping overall frequencies",
    addFreqOverall = "Calculating overall frequencies",

    # Imports
    singleVcf = "Importing from single VCF.",
    multiVcfs = "Importing from multiple VCF files."
)

# EnsDb packages ----------------------------------------------------------

## list all packages...
EnsDbPacks <- grep(
    "^EnsDb", rownames(installed.packages()), value = TRUE)
