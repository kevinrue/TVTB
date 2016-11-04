
stopifnot(
    require(TVTB),
    require(BiocParallel),
    require(rtracklayer),
    require(Rsamtools),
    require(VariantAnnotation),
    require(ensemblVEP),
    require(ggplot2),
    require(DT),
    require(reshape2),
    require(dplyr)
)

message(
    "Please run ",
    "devtools::install_github(c(\"rstudio/shiny\", \"rstudio/DT\")) ",
    "to enable latest functionalities.")

# Seems like all the following don't need to be re-imported for the app
# library(GenomicRanges)
# library(reshape2)
# library(ensembldb)

# General settings --------------------------------------------------------

# Alphabetical order (ignoring ^default, ^choices, ^all, ...)
GS <- list(
    default.vep = "CSQ",
    choices.vepCols = c(),
    default.vepCols = c(),
    default.refGenotypes = c("0|0"),
    default.hetGenotypes = c("0|1", "1|0", "1|2", "2|1"),
    default.altGenotypes = c("1|1", "2|2"),
    choices.ensDbType = list(
        "Gene name" = "Genename"),
    default.ensDbType = "gene",
    choices.ensDbFilters = c("=", "!=", "like", "in"),
    default.ensDbFilters = "=",
    all.genotypes = c(
        "0/0","0/1","1/0","1/1","./0","0/.","./1","1/.","./.",
        "0/2","2/0", "1/2", "2/1", "2/2", "./2", "2/.",
        "0|0","0|1","1|0","1|1",".|0","0|.",".|1","1|.",
        "0|2","2|0", "1|2", "2|1", "2|2", ".|2", "2|.",
        "0", "1", "2","."),
    choices.phenoCols = c(),
    default.phenoCols = c(),
    choices.regionInputMode = list(
        "BED file" = "bed",
        "UCSC browser" = "ucsc",
        "EnsDb package" = "EnsDb"
    ),
    default.regionInputMode = c("BedFile"),
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

# Display settings ----

vepCountHeight <- "500px"

# R session variables -----------------------------------------------------

originalOptions <- options()
options("width" = 120)

# Messages ----------------------------------------------------------------

# Alphabetical order (ignoring ^default, ^choices, ^all, ...)
Msgs <- list(
    # Mandatory inputs
    annotationPackage = "An EnsDb annotation package must be selected.",
    hetGenotypes = "Heterozygote genotype(s) must be defined.",
    altGenotypes = "Alternative homozygote genotype(s) must be defined.",
    vepKey = "INFO field of VEP prediction must be defined.",
    maf.min = "Maximum minor allele frequency",
    maf.min = "Mimimum minor allele frequency",
    refGenotypes = "Reference genotype(s) must be defined",
    invalidUcscRegions = "Invalid UCSC-type input.",
    vcfFolder = "VCF folder",
    vcfPattern = "VCF pattern",
    xAxisAngle = "X axis angle",
    xAxisHjust = "X axis horizontal justification",
    xAxisVjust = "X axis vertical justification",
    xAxisSize = "X axis font size",

    # Tabulate VEP by phenotype (TVBP)
    vepTVBP = "VEP field tabulated",
    vepFacetKeyTVBP = "Faceting VEP field",
    stackedPercentageTVBP = "Stack percentage?",

    # Density VEP by phenotype (DVBP)
    vepDVBP = "VEP field tabulated",
    vepFacetKeyDVBP = "Faceting VEP field",

    # Optional inputs
    importVariants = "Please import/refresh variants.",
    genomicRanges = "No genomic range defined. All variants considered.",
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
    fileChooseCancelled = "file choice cancelled."
)

# Tracking messages -------------------------------------------------------

# Progress bars
Tracking = list(
    calculate = "Crunching data...",
    ggplot = "Assembling plot (ggplot2)",
    render = "Rendering.",
    preprocessing = "Pre-processing data",
    postprocessing = "Post-processing data",
    singleVcf = "Importing from single VCF.",
    multiVcfs = "Importing from multiple VCF files."
)

# EnsDb packages ----------------------------------------------------------

## list all packages...
packs <- rownames(installed.packages())
EnsDbPacks <- packs[grep(packs, pattern="^EnsDb")]
