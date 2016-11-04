\name{vepInPhenoLevel-methods}

\alias{vepInPhenoLevel}
\alias{vepInPhenoLevel,ExpandedVCF,TVTBparam-method}
\alias{vepInPhenoLevel,ExpandedVCF,missing-method}

\docType{methods}
\alias{vepInPhenoLevel-methods}
\title{VEP predictions of variants observed in samples}
\description{
Considers only variants observed (uniquey) in samples associated with a given
phenotype level, and tabulates the corresponding values for a given VEP
prediction field.
}

\usage{
\S4method{vepInPhenoLevel}{ExpandedVCF,TVTBparam}(
    vcf, phenoCol, level, vepCol, param, ...,
    unique = FALSE, facet = NULL)

\S4method{vepInPhenoLevel}{ExpandedVCF,missing}(
    vcf, phenoCol, level, vepCol, alts, param = NULL, ...,
    unique = FALSE, facet = NULL)
}

\arguments{
    \item{vcf}{
\code{ExpandedVCF} object.
}
    \item{phenoCol}{
Name of a column in \code{pheno}.
}
    \item{level}{
Phenotype level; only variants observed in at least one sample will be
considered.
}
    \item{vepCol}{
VEP prediction field; Name of a metadata column in
\code{ensemblVEP::parseCSQToGRanges(vcf)}.
}
    \item{param}{
\code{TVTBparam} object that contains recurrent parameters.
}
    \item{alts}{
\code{character} vector of alternate genotypes
(both heterozygote and homozygote).
}
    \item{...}{
Additional arguments used to override \code{TVTBparam} slots.
}
    \item{unique}{
If \code{TRUE}, consider only variants unique to the phenotype level
(\emph{i.e.} absent from all other phenotype levels).
}
    \item{facet}{
Name of a metadata column in \code{ensemblVEP::parseCSQToGRanges(vcf)}.
Additional VEP field appended as an additional column to the
\code{data.frame} returned. If \code{plot=TRUE}, this field will be used to
create one sub-plot for each level of the faceting field.

If available, \code{"Feature"} is a recommended value for this argument,
as VEP typically produce one prediction per variant per feature.
}
}

\value{
A \code{data.frame} in long format including one row per variant
seen in at least one sample (heterozygote or alternate homozygote)
associated with the phenotype level.
The \code{data.frame} contains at least one column for the VEP prediction
value.
An additional column containing another VEP prediction field may be added
using the \code{facet} argument.
}
\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\linkS4class{VCF}},
\code{\link{ensemblVEP}},
\code{\linkS4class{GRanges}},
and \code{\linkS4class{DataFrame}}.
}

\examples{
# Example data ----

# VCF file
extdata <- file.path(system.file(package = "TVTB"), "extdata")
vcfFile <- file.path(extdata, "moderate.vcf")

# Phenotype file
phenoFile <- file.path(extdata, "moderate_pheno.txt")
phenotypes <- S4Vectors::DataFrame(
    read.table(file = phenoFile, header = TRUE, row.names = 1))

# TVTB parameters
tparam <- TVTBparam(
    genos = list(
        REF = "0|0",
        HET = c("0|1", "1|0"),
        ALT = "1|1"))

# Pre-process variants ----

# Import variants
vcf <- VariantAnnotation::readVcf(file = vcfFile)
# Add phenotype information necessary for the demo
colData(vcf) <- phenotypes
# Separate multi-allelic records into bi-allelic records
vcf <- VariantAnnotation::expand(x = vcf, row.names = TRUE)
# Disambiguate row.names from multi-allelic records
rownames(vcf) <- paste(rownames(vcf), mcols(vcf)[,"ALT"], sep = "_")

# Example usage ----

vepInPhenoLevel(
    vcf = vcf, phenoCol = "super_pop", level = "AFR",
    vepCol = "CADD_PHRED", param = tparam)
}

\keyword{methods}

