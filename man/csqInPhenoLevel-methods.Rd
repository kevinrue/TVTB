\name{csqInPhenoLevel-methods}

\alias{csqInPhenoLevel}
\alias{csqInPhenoLevel,ExpandedVCF,TVTBparam-method}
\alias{csqInPhenoLevel,ExpandedVCF,missing-method}

\docType{methods}
\alias{csqInPhenoLevel-methods}
\title{VEP predictions of variants observed in samples}
\description{
Considers only variants observed (uniquey) in samples associated with a given
phenotype level, and tabulates the corresponding values for a given VEP
prediction field.
}

\usage{
\S4method{csqInPhenoLevel}{ExpandedVCF,TVTBparam}(
    vcf, phenoCol, level, csqCol, param, ...,
    unique = FALSE, facet = NULL)

\S4method{csqInPhenoLevel}{ExpandedVCF,missing}(
    vcf, phenoCol, level, csqCol, alts, param = NULL, ...,
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
    \item{csqCol}{
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

# Parameters
tparam <- TVTBparam(
    genos = list(
        REF = "0|0",
        HET = c("0|1", "1|0"),
        ALT = "1|1"))

vcf <- preprocessVariants(
    file = vcfFile, param = tparam, phenos = phenoFile)

# Example usage ----

csqInPhenoLevel(
    vcf = vcf, phenoCol = "super_pop", level = "AFR",
    csqCol = "CADD_PHRED", param = tparam)
}

\keyword{methods}

