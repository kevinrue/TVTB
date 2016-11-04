library(TVTB)

## ----TVTBparamCreate---------------------------------------------
tparam <- TVTBparam(
    genos = list(
        REF = c("0|0"),
        HET = c("0|1", "1|0"),
        ALT = c("1|1")),
    ranges = GenomicRanges::GRanges(
        seqnames = "15",
        IRanges::IRanges(
            start = 48413170, end = 48434757, names = "SLC24A5")),
    aaf = "AAF", # default
    maf = "MAF", # default
    vep = "CSQ", # default
    bp = BiocParallel::SerialParam()) # default

## ----TVTBparamView-----------------------------------------------
tparam

## ----importPhenotypes--------------------------------------------
extdata <- system.file("extdata", package = "TVTB")
phenoFile <- file.path(extdata, "integrated_samples.txt")
phenotypes <- S4Vectors::DataFrame(
    read.table(file = phenoFile, header = TRUE, row.names = 1))

## ----vcfFile-----------------------------------------------------
vcfFile <- file.path(extdata, "chr15.phase3_integrated.vcf.gz")
tabixVcf <- Rsamtools::TabixFile(file = vcfFile)

## ----ScanVcfParam------------------------------------------------
svp <- VariantAnnotation::ScanVcfParam(
    fixed = c("ALT", "QUAL", "FILTER"), # all fields: could be omitted
    info = c(vep(tparam)),
    geno = "GT", # all fields: could be omitted
    samples = rownames(phenotypes),
    which = ranges(tparam))

## ----preprocessVariants, message=FALSE---------------------------
library(SummarizedExperiment)
# Import variants as a CollapsedVCF object
vcf <- VariantAnnotation::readVcf(file = vcfFile, param = svp)
# Combine with phenotype information
colData(vcf) <- phenotypes
# Expand into a ExpandedVCF object (bi-allelic records)
vcf <- VariantAnnotation::expand(x = vcf, row.names = TRUE)

## ----addOverallFrequencies, message=FALSE------------------------
vcf <- addOverallFrequencies(vcf = vcf, param = tparam)
vcf <- addPhenoLevelFrequencies(
    vcf = vcf, pheno = "super_pop", level = "AFR", param = tparam)

## ----addOverallFrequencies, message=FALSE------------------------
fixedR <- VcfFixedRules(exprs = list(
    pass = expression(FILTER == "PASS"),
    qual = expression(QUAL > 20)
))
fixedR

infoR <- VcfInfoRules(exprs = list(
    common = expression(MAF > 0.1), # minor allele frequency
    present = expression(ALT + HET > 0) # count of non-REF homozygotes
))
# ...is synonym to...
infoR <- VcfInfoRules(exprs = list(
    common = expression(MAF > 0.1), # minor allele frequency
    present = expression(ALT > 0 | HET > 0)
))
infoR

vepR <- VcfVepRules(exprs = list(
    missense = expression(Consequence %in% c("missense_variant")),
    CADD = expression(CADD_PHRED > 15)
))
vepR

vcfRules <- VcfFilterRules(fixedR, infoR, vepR)
vcfRules

active(vcfRules)["CADD"] <- FALSE

eval(expr = vcfRules, envir = vcf)
