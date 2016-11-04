# This script was used to:
## import variants from the 1000 genomes,
## calculate and add genotype frequencies to the VCF object,
## subset variants to only MODERATE IMPACT VEP predictions
## subset samples to only "MSL", "GBR" populations

# Links:
## ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

# Import and clean data ----

# Package data folder
extdata <- file.path(system.file(package = "TVTB"), "extdata")

# VCF file (1000 genomes)
vcfFile <- file.path(extdata, "chr15.phase3_integrated.vcf.gz")
tabixVcf <- Rsamtools::TabixFile(file = vcfFile)

# Phenotypes (1000 genomes)
phenoFile <- file.path(extdata, "integrated_samples.txt")
phenotypes <- S4Vectors::DataFrame(read.table(
    file = phenoFile, header = TRUE, row.names = 1))

# VCF params
svp <- VariantAnnotation::ScanVcfParam(
    info = "CSQ",
    samples = rownames(phenotypes))

# Import variants in a CollapsedVCF object
vcf <- VariantAnnotation::readVcf(file = tabixVcf, param = svp)

# Attach the phenotype data to the CollapsedVCF object
colData(vcf) <- phenotypes

# Trim from the header INFO fields without imported data
info(vcf) <- subset(info(vcf), select = colnames(info(vcf)) == "CSQ")
info(header(vcf)) <- subset(
    info(header(vcf)), subset = rownames(info(header(vcf))) == "CSQ")

# Extract Ensembl VEP predictions
csq <- ensemblVEP::parseCSQToGRanges(x = vcf, VCFRowID = rownames(vcf))

# Reduce number of variants ----

# How many of each IMPACT predictions?
table(mcols(csq)[,"IMPACT"])

# Extract the variants associated with the 62 MODERATE impact predictions
vcf.moderate <- vcf[
    unique(subset(
        x = mcols(csq),
        subset = mcols(csq)[,"IMPACT"] == "MODERATE",
        select = "VCFRowID",
        drop = TRUE))
    ,]

# Reduce number of samples ----

# Extract genotype matrix for MODERATE variants
GT <- geno(vcf.moderate)[["GT"]]
# List of genotypes observed in the matrix
genoNames <- unique(as.vector(GT))
# Define carrier genotypes
altGenos <- grep(
    pattern = "0|0", x = genoNames, invert = TRUE, fixed = TRUE, value = TRUE)

## How many variants each sample carries ?
# apply(X = GT, MARGIN = 2, FUN = function(x)(sum(x %in% altGenos)))
## How many samples carry each count of variants ?
# table(apply(X = GT, MARGIN = 2, FUN = function(x)(sum(x %in% altGenos))))

# Function to return whether a variant is observed in each level of
# a phenotype
carrierByPhenoLevel <- function(vcf, rowIdx, pheno, altGenos){
    tapply(
        X = geno(vcf)[["GT"]][rowIdx,],
        INDEX = colData(vcf)[,pheno],
        FUN = function(x)(sum(x %in% altGenos) > 0))
}

# Function to return whether each variant is observed in each level of
# of a phenotype
carrierByPhenoLevelByVariant <- function(vcf, pheno, altGenos){
    t(sapply(
        X = seq_along(vcf),
        FUN = function(x)(
            carrierByPhenoLevel(
                vcf = vcf, rowIdx = x, pheno = pheno, altGenos = altGenos)
            )
        ))
}

# Which variants are seen in each population?
carriers.pop <- carrierByPhenoLevelByVariant(
    vcf = vcf.moderate, pheno = "pop", altGenos = altGenos)
# Total carrier by pop, by super_pop
sapply(
    X = levels(phenotypes$super_pop),
    FUN = function(x){
        pop <- unique(subset(phenotypes, subset = super_pop == x, select = "pop", drop = TRUE))
        colSums(carriers.pop)[pop]
        })

# Which variants are seen in each super-population?
carriers.super_pop <- carrierByPhenoLevelByVariant(
    vcf = vcf.moderate, pheno = "super_pop", altGenos = altGenos)
colSums(carriers.super_pop)

# Keep only AFR_MSL (8 carriers) and EUR_GBR (0 carriers) as examples
vcf.moderate.pop <- vcf.moderate[,which(phenotypes$pop %in% c("MSL", "GBR"))]
# Write the subsetted VCF file
writeVcf(
    obj = vcf.moderate.pop,
    filename = "~/Dropbox/TVTB/inst/extdata/moderate.vcf")

# Keep only phenotype information for the relevant samples
phenotypes.pop <- subset(phenotypes, pop %in% c("MSL", "GBR"))
# Write the subsetted phenotype information
write.table(
    x = cbind(
        sample = rownames(phenotypes.pop),
        as.data.frame(phenotypes.pop)),
    file = "~/Dropbox/TVTB/inst/extdata/moderate_pheno.txt",
    quote = FALSE,
    row.names = FALSE)
