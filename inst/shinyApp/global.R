
stopifnot(
    require(tSVE),
    require(BiocParallel),
    require(rtracklayer),
    require(Rsamtools),
    require(VariantAnnotation),
    require(ensemblVEP),
    require(ggplot2),
    require(DT)
)

message(
    "Please make sure 'DT' package was installed using ",
    "devtools::install_github('rstudio/DT') ",
    "to enable latest functionalities.")

# Seems like all the following don't need to be re-imported for the app
# library(GenomicRanges)
# library(reshape2)
# library(ensembldb)

# General settings --------------------------------------------------------

# Alphabetical order (ignoring ^default, ^choices, ^all, ...)
GS <- list(
    default.csq = "CSQ",
    choices.csqCols = c(),
    default.csqCols = c(),
    default.refGenotypes = c("0|0"),
    default.hetGenotypes = c("0|1", "1|0"),
    default.altGenotypes = c("1|1"),
    choices.ensDbType = list(
        "Gene name" = "Genename"),
    default.ensDbType = "gene",
    choices.ensDbFilters = c("=", "!=", "like", "in"),
    default.ensDbFilters = "=",
    choices.genome = list("Auto (VCF header)" = "auto", "hg19"),
    default.genome = "auto",
    all.genotypes = c(
        "0/0","0/1","1/0","1/1","./0","0/.","./1","1/.","./.",
        "0|0","0|1","1|0","1|1",".|0","0|.",".|1","1|.",
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
    default.vcfFolder = system.file("extdata", package = "tSVE"),
    default.vcfPattern = "^chr%s\\..*\\.vcf\\.gz$",
    choices.vcfInputMode = list(
        "Single VCF" = "SingleVcf",
        "One per chromosome" = "OnePerChr"),
    default.vcfInputMode = "OnePerChr"
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

csqCountHeight <- "500px"

# Messages ----------------------------------------------------------------

# Alphabetical order (ignoring ^default, ^choices, ^all, ...)
Msgs <- list(
    # Mandatory inputs
    annotationPackage = "Annotation package",
    hetGenotypes = "Heterozygote genotype(s)",
    altGenotypes = "Alternative homozygote genotype(s)",
    bedFile = "BED file",
    bedRecords = "Genomic regions must be provided",
    bpConfig = "BiocParallel configuration",
    csqCols = "Fields of VEP predictions",
    csqField = "VEP INFO field",
    doGenoHeatmap = "Please click the button to generate/update the figure",
    ensDb.type = "Type of the Ensembl filter",
    ensDb.condition = "Condition of the Ensembl filter",
    ensDb.values = "Values for the Ensembl filter",
    facet = "Faceting VEP field",
    genoFirstCol = "Index of first column in genotype matrix",
    genoFirstRow = "Index of first row of genotype matrix",
    genoNumCols = "Number of columns in genotype matrix",
    genoNumRows = "Number of rows in genotype matrix",
    maf.min = "Maximum minor allele frequency",
    maf.min = "Mimimum minor allele frequency",
    phenoCols = "Phenotype fields",
    phenotype = "Phenotype field",
    refGenotypes = "Reference genotype(s)",
    refreshVariants = "Please refresh variants",
    selectBed = "Please select BED file using the 'Browse' button.",
    selectVcf = "Please select VCF file using the 'Browse' button.",
    stackedPercentage = "Stack percentage?",
    ucscRegions = "UCSC-type region(s)",
    variantCsq = "Variant prediction field",
    vcfCols = "VCF columns",
    vcfFolder = "VCF folder",
    vcfPattern = "VCF pattern",
    xAxisAngle = "X axis angle",
    xAxisHjust = "X axis horizontal justification",
    xAxisVjust = "X axis vertical justification",
    xAxisSize = "X axis font size",

    # Optional inputs
    phenoFile = "A Phenotype file may be provided",

    # Catch warnings/errors
    vcf = "Impossible to read VCF file. See console.",
    vcfs = "Impossible to read VCF file. See console.",

    # Calculated values (missing are likely bugs)
    bedChrs = "bedChrs bug?",
    bpParam = "bpParam bug?",
    chrVcf = "chrVcf bug?",
    csq = "csq bug?",
    csqMcols = "csqMcols bug?",
    csqTable = "csqTable bug?",
    csqTableTable = "csqTableTable bug?",
    csqTableTableDecreasing = "csqTableTableDecreasing bug?",
    edb = "edb bug?",
    ensDbFilter = "Invalid EnsDb filter",
    genomeSeqinfo = "genomeSeqinfo bug?",
    genoSampleRanges.rows = "genoSampleRanges$rows bug?",
    genoSampleRanges.cols = "genoSampleRanges$cols bug?",
    genotypes = "genotypes bug?",
    legendText = "legendText bug?",
    mafRange = "mafRange bug?",
    phenotypes = "phenotypes bug?",
    phenos = "phenos bug?",
    queryGenes = "queryGenes bug?",
    singleVcf = "singleVcf bug?",
    tfl = "tfl bug?",
    tParam = "tParam bug?",
    varCsqPlotX = "varCsqPlotX bug?",
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
    checkPhenotypes = "Checking VCF samples against phenotype file",
    singleVcf = "Importing from single VCF.",
    multiPaths = "Checking paths to multiple VCFs",
    multiVcfs = "Importing from multiple VCF files.",
    mergeVcfs = "Merging multiple VCF objects"
)

# EnsDb packages ----------------------------------------------------------

## list all packages...
packs <- installed.packages()
EnsDbPacks <- packs[grep(packs, pattern="^EnsDb")]

# R session variables -----------------------------------------------------

options("width" = 120)
