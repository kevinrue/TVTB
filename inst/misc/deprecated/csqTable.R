require(EnsDb.Hsapiens.v75)
library(GenomicRanges)
library(IRanges)
library(VariantAnnotation)
library(Rsamtools)
library(GenomeInfoDb)
library(ensemblVEP)
library(BiocParallel)

# Example data ----
extdata <- file.path(system.file(package = "tSVE"), "extdata")
vcfFile <- file.path(extdata, "chr1_sample.vcf.gz")
phenoFile <- file.path(extdata, "pheno_sample.txt")
bedRegions <- GRanges(
    seqnames = 1,
    ranges = IRanges(start = 211744929, end = 211745433))

# Import phenotype data ----
phenotypes <- read.table(
    file = phenoFile,
    header = TRUE,
    row.names = 1)

# Import variants ----
svp <- ScanVcfParam(
    geno = "GT",
    samples = rownames(phenotypes),
    which = bedRegions)
tf <- TabixFile(vcfFile)
si <- seqinfo(EnsDb.Hsapiens.v75)
vcf <- readVcf(file = tf, genome = si, param = svp)
csq <- parseCSQToGRanges(
    vcf,
    VCFRowID = rownames(vcf),
    info.key = "ANN")

# Example usage ----
tabulateCsqPhenotype(
    level = "0", phenos = as.factor(phenotypes[,"PAH"]),
    genos = geno(vcf)[["GT"]],
    altGenos = c("0/1", "1/1"), csqs = csq, csq.field = "Consequence",
    unique = FALSE, facet = NULL)

tabulateCsqPhenotypes(
    phenos = phenotypes[,"PAH"], genos = geno(vcf)[["GT"]],
    altGenos = c("0/1", "1/1"), csqs = csq, csq.field = "Consequence",
    unique = FALSE, facet = NULL, BPPARAM = SerialParam())

csqTable <- tabulateCsqPhenotypes(
    phenos = phenotypes[,"PAH"], genos = geno(vcf)[["GT"]],
    altGenos = c("0/1", "1/1"), csqs = csq, csq.field = "Consequence",
    unique = FALSE, facet = "Feature", BPPARAM = SerialParam())

csqTableTable <- as.data.frame(table(csqTable))

ggplot(data = csqTable, mapping = aes(Phenotype, fill = "Consequence")) +
    scale_x_discrete(drop = FALSE) +
    geom_bar() +
    facet_wrap(facets = "Feature")
