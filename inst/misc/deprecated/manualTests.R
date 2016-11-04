
library(IRanges)
library(GenomicRanges)
library(Rsamtools)
library(VariantAnnotation)
library(EnsDb.Hsapiens.v75)
library(ensemblVEP)
library(BiocParallel)

extdata <- "inst/extdata/"
vcfFile <- "chr1_sample.vcf.gz"
phenoFile <- "pheno_sample.txt"
bedRegions <- GRanges(
    seqnames = 1,
    ranges = IRanges(start = 211744929, end = 211745433))

phenotypes <- read.table(
    file = file.path(extdata, phenoFile),
    header = TRUE,
    row.names = 1)

svp <- ScanVcfParam(
    geno = "GT",
    samples = rownames(phenotypes),
    which = bedRegions)
tf <- TabixFile(file.path(extdata, vcfFile))
si <- seqinfo(EnsDb.Hsapiens.v75)
vcf <- readVcf(file = tf, genome = si, param = svp)
csq <- parseCSQToGRanges(
    vcf,
    VCFRowID = rownames(vcf),
    info.key = "ANN")


tabulatePredictionsPhenotype(
    l = "0", phenos = as.factor(phenotypes[,"PAH"]), genos = geno(vcf)[["GT"]],
    altGenos = c("0/1", "1/1"), csqs = csq, csq.field = "Consequence",
    unique = FALSE, facet = NULL)

csqTableByPheno(
    phenos = phenotypes[,"PAH"], genos = geno(vcf)[["GT"]],
    altGenos = c("0/1", "1/1"), csqs = csq, csq.field = "Consequence",
    unique = FALSE, facet = NULL, BPPARAM = SerialParam())
