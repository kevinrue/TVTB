
library(VariantAnnotation)

# Example data ----

extdata <- file.path(system.file(package = "tSVE"), "extdata")
vcfFile <- file.path(extdata, "chr15.phase3_integrated.vcf.gz")
phenoFile <- file.path(extdata, "integrated_samples.txt")
bedRegions <- GenomicRanges::GRanges(
    seqnames = "15",
    ranges = IRanges::IRanges(start = 48413169, end = 48434869))

# Import phenotype data ----

phenotypes <- S4Vectors::DataFrame(read.table(
    file = phenoFile,
    header = TRUE,
    row.names = 1))

# Import variants ----

svp <- VariantAnnotation::ScanVcfParam(
    geno = "GT",
    info = "CSQ",
    samples = rownames(phenotypes),
    which = bedRegions)
tf <- Rsamtools::TabixFile(vcfFile)
#
vcf <- VariantAnnotation::readVcf(file = tf, param = svp)
# Expand pontentially multi-allelic records to bi-allelic records
eVcf <- VariantAnnotation::expand(x = vcf, row.names = TRUE)
# Disambiguate row.names of multi-allelic records
rownames(eVcf) <- paste(rownames(eVcf), mcols(eVcf)[,"ALT"], sep = "_")

# Example usage ----
variantsInSamples(
    vcf = eVcf,
    altGenos = c("0|1", "1|0", "1|1"),
    samplesIdx = which(phenotypes$super_pop == "EUR"),
    unique = FALSE)

# Subset multi-allelic records
vcfMulti <- vcf[which(lengths(mcols(vcf)[,"ALT"]) > 1)]

# Which genotypes are present?
table(geno(vcfMulti)[["GT"]][1,])

#  0|0  0|1  0|2  1|0  2|0
# 2496    1    3    2    2

# As a consequence, users may ask "which variants are present in any ALT form"
