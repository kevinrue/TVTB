\name{vepInPhenoLevel-methods}

\alias{vepInPhenoLevel}
\alias{vepInPhenoLevel,ExpandedVCF-method}

\docType{methods}
\alias{vepInPhenoLevel-methods}
\title{VEP predictions of variants observed in samples}
\description{
Returns VEP predictions for variants observed (uniquey) in samples associated
with a given phenotype level.
}

\usage{
\S4method{vepInPhenoLevel}{ExpandedVCF}(
    vcf, phenoCol, level, vepCol, unique = FALSE)
}

\arguments{
    \item{vcf}{
\code{ExpandedVCF} object.

\code{metadata(vcf)[["TVTBparam"]]} must contain a
\code{\linkS4class{TVTBparam}} object.
}
    \item{phenoCol}{
Name of a column in \code{pheno}.
}
    \item{level}{
Phenotype level; only variants observed in at least one sample will be
considered.
}
    \item{vepCol}{
VEP prediction fields; \code{character} vector of metadata columns in
\code{ensemblVEP::parseCSQToGRanges(vcf)}.
}
    \item{unique}{
If \code{TRUE}, consider only variants unique to the phenotype level
(\emph{i.e.} absent from all other phenotype levels).
}
}

\value{
A \code{GRanges} including all VEP predictions associated with a variant
seen in at least one sample (heterozygote or alternate homozygote)
associated with the phenotype level.
The \code{GRanges} contains at least one column for the VEP prediction
value.
Additional columns containing another VEP prediction field may be added
using the \code{facet} argument.
}

\section{Note}{
If available, \code{"Feature"} is a recommended value for this argument,
as VEP typically produce one prediction per variant per feature.
}

\section{Warning}{
A \code{warning} message is issued if genotypes are not fully defined in the
\code{TVTBparam}.
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
vcfFile <- system.file("extdata", "moderate.vcf", package = "TVTB")

# Phenotype file
phenoFile <- system.file("extdata", "moderate_pheno.txt", package = "TVTB")
phenotypes <- S4Vectors::DataFrame(
    read.table(file = phenoFile, header = TRUE, row.names = 1))

# TVTB parameters
tparam <- TVTBparam(Genotypes("0|0", c("0|1", "1|0"), "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(
    vcfFile, param = tparam, colData = phenotypes)
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)


# Example usage ----

vepInPhenoLevel(vcf, "super_pop", "AFR", c("CADD_PHRED", "Feature", "IMPACT"))
}

\keyword{methods}
