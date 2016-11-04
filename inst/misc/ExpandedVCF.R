
extdata <- system.file("extdata", package = "TVTB")
list.files(extdata)
vcfFile <- file.path(extdata, "chr15.phase3_integrated.vcf.gz")

library(S4Vectors)
phenoFile <- file.path(extdata, "integrated_samples.txt")
phenotypes <- read.table(file = phenoFile, header = TRUE, row.names = 1)

library(VariantAnnotation)
vcf <- readVcf(file = vcfFile)


colData(vcf) <- DataFrame(phenotypes)
