extdata <- file.path(system.file(package = "tSVE"), "extdata")

# VCF file
vcfFile <- file.path(extdata, "chr15.phase3_integrated.vcf.gz")
tabixVcf <- Rsamtools::TabixFile(file = vcfFile)

# Pheno
phenoFile <- file.path(extdata, "integrated_samples.txt")
phenotypes <- S4Vectors::DataFrame(read.table(
    file = phenoFile, header = TRUE, row.names = 1))

# VCF params
svp <- VariantAnnotation::ScanVcfParam(
    info = "CSQ",
    samples = rownames(phenotypes))

# VCF
vcf <- VariantAnnotation::readVcf(file = tabixVcf, param = svp)

# Combine the phenotype data
colData(vcf) <- phenotypes

# Trim unused INFO fields
info(vcf) <- subset(info(vcf), select = colnames(info(vcf)) == "CSQ")
info(header(vcf)) <- subset(info(header(vcf)), subset = rownames(info(header(vcf))) == "CSQ")

# Extract consequences
csq <- ensemblVEP::parseCSQToGRanges(x = vcf, VCFRowID = rownames(vcf))

# How many of each IMPACT predictions? 62
table(mcols(csq)[,"IMPACT"])

# Extract the variants associated with the 62 MODERATE impact predictions
vcf.moderate <- vcf[
    unique(subset(
        x = mcols(csq),
        subset = mcols(csq)[,"IMPACT"] == "MODERATE",
        select = "VCFRowID",
        drop = TRUE))
    ,]

# How many samples of each pop phenotype carry those variants?
GT <- geno(vcf.moderate)[["GT"]]
genoNames <- unique(as.vector(GT))
altGenos <- grep(pattern = "0|0", x = genoNames, invert = TRUE, fixed = TRUE, value = TRUE)

# apply(X = GT, MARGIN = 2, FUN = function(x)(sum(x %in% altGenos)))
# table(apply(X = GT, MARGIN = 2, FUN = function(x)(sum(x %in% altGenos))))

carrierByPhenoLevel <- function(vcf, rowIdx, pheno, altGenos){
    tapply(
        X = geno(vcf)[["GT"]][rowIdx,],
        INDEX = colData(vcf)[,pheno],
        FUN = function(x)(sum(x %in% altGenos) > 0))
}

carrierByPhenoLevelByVariant <- function(vcf, pheno, altGenos){
    t(sapply(
        X = seq_along(vcf),
        FUN = function(x)(
            carrierByPhenoLevel(
                vcf = vcf, rowIdx = x, pheno = pheno, altGenos = altGenos)
            )
        ))
}

carriers.pop <- carrierByPhenoLevelByVariant(
    vcf = vcf.moderate, pheno = "pop", altGenos = altGenos)
# Total carrier by pop, by super_pop
sapply(
    X = levels(phenotypes$super_pop),
    FUN = function(x){
        pop <- unique(subset(phenotypes, subset = super_pop == x, select = "pop", drop = TRUE))
        colSums(carriers.pop)[pop]
        })

carriers.super_pop <- carrierByPhenoLevelByVariant(
    vcf = vcf.moderate, pheno = "super_pop", altGenos = altGenos)
colSums(carriers.super_pop)

# Keep only AFR_MSL (8 carriers) and EUR_GBR (0 carriers)
vcf.moderate.pop <- vcf.moderate[,which(phenotypes$pop %in% c("MSL", "GBR"))]

writeVcf(
    obj = vcf.moderate.pop,
    filename = "~/Dropbox/tSVE/inst/extdata/moderate.vcf")

phenotypes.pop <- subset(phenotypes, pop %in% c("MSL", "GBR"))
write.table(
    x = cbind(
        sample = rownames(phenotypes.pop),
        as.data.frame(phenotypes.pop)),
    file = "~/Dropbox/tSVE/inst/extdata/moderate_pheno.txt",
    quote = FALSE,
    row.names = FALSE)
