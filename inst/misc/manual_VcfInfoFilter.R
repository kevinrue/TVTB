library(VariantAnnotation)

# Example data ----

extdata <- file.path(system.file(package = "tSVE"), "extdata")
vcfFile <- file.path(extdata, "moderate.vcf")
phenoFile <- file.path(extdata, "moderate_pheno.txt")

# Import phenotype data ----

phenotypes <- S4Vectors::DataFrame(read.table(
    file = phenoFile,
    header = TRUE,
    row.names = 1))

# Define tSVEParam

tparam <- tSVEParam(
    genos = list(
        REF = "0|0",
        HET = c("0|1", "1|0"),
        ALT = "1|1"))

# Example usage ----

preprocessVariants(
    file = vcfFile, param = tparam, phenos = phenotypes)

# Add overall frequencies
vcf <- addOverallFrequencies(vcf = vcf, param = tparam)

#
vif <- VcfInfoFilter(name = "MAF", condition = ">", value = 0.01)

# Only one variant MAF > 0

vcfTest <- parse(text = sprintf(
    "info(vcf)[,\"%s\"] %s %s",
    slot(vif, "name"),
    slot(vif, "condition"),
    slot(vif, "value")))

eval(vcfTest)

