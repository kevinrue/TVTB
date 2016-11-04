\name{densityCsq-methods}

\alias{densityCsq}

\alias{densityCsqByPhenotype}
\alias{densityCsqInPhenoLevel}

\alias{densityCsqByPhenotype,ExpandedVCF,tSVEParam-method}
\alias{densityCsqByPhenotype,ExpandedVCF,missing-method}

\alias{densityCsqInPhenoLevel,ExpandedVCF,tSVEParam-method}
\alias{densityCsqInPhenoLevel,ExpandedVCF,missing-method}

\docType{methods}
\alias{densityCsq-methods}
\alias{densityCsqByPhenotype-methods}
\alias{densityCsqInPhenoLevel-methods}
\title{Density of VEP predictions in phenotypes}
\description{
\code{densityCsqByPhenotype} plots the density of a given VEP prediction
    field (coerced to \code{numeric}) for variants observed (uniquely) in each
    level of a given phenotype.

\code{densityCsqInPhenoLevel} considers only variants observed in samples
    associated with a given phenotype level, and plots the corresponding
    density for a given VEP prediction.
}

\usage{
\S4method{densityCsqByPhenotype}{ExpandedVCF,tSVEParam}(
        vcf, phenoCol, csqCol, param, ...,
        unique = FALSE, facet = NULL, plot = FALSE, popFreq = FALSE)

\S4method{densityCsqByPhenotype}{ExpandedVCF,missing}(
        vcf, phenoCol, csqCol, alts, param = NULL, ...,
        unique = FALSE, facet = NULL, plot = FALSE, popFreq = FALSE)

\S4method{densityCsqInPhenoLevel}{ExpandedVCF,tSVEParam}(
        level, vcf, phenoCol, csqCol, param, ...,
        unique = FALSE, facet = NULL, plot = FALSE, popFreq = FALSE)

\S4method{densityCsqInPhenoLevel}{ExpandedVCF,missing}(
        level, vcf, phenoCol, csqCol, alts, param = NULL, ...,
        unique = FALSE, facet = NULL, plot = FALSE, popFreq = FALSE)
}

\arguments{
    \item{vcf}{
\code{ExpandedVCF} object.
}
    \item{phenoCol}{
Name of a column in \code{pheno}.
}
    \item{csqCol}{
Name of a metadata column in \code{ensemblVEP::parseCSQToGRanges(vcf)}.
}
    \item{param}{
\code{tSVEParam} object that contains recurrent parameters.
}
    \item{alts}{
Character vector of alternate genotypes (both heterozygote and homozygote).
}
    \item{level}{
Phenotype level; only variants observed in at least one sample will be
considered.
}
    \item{...}{
Additional arguments used to override \code{tSVEParam} slots.
}
    \item{unique}{
If \code{TRUE}, consider only variants unique to the phenotype level
(\emph{i.e.} not seen in any other phenotype level).
}
    \item{facet}{
Name of a column in \code{ensemblVEP::parseCSQToGRanges(vcf)}.
Additional VEP field appended as an additional column to the
\code{data.frame} returned. If \code{plot=TRUE}, this field will be used to
create one sub-plot for each level of the faceting field.

If available, \code{"Feature"} is a recommended value for this argument,
as VEP typically produce one prediction per variant per feature.
}
    \item{plot}{
If \code{TRUE}, return a \code{ggplot} object.
}
    \item{popFreq}{
Should be set to \code{TRUE} if \code{csqField} is a VEP field representing a
population frequency of the form \code{allele:frequency}.
If set to \code{TRUE}, only \code{frequency} will be extracted and typed as
\code{numeric}.
}
}

\value{
\describe{
\item{\code{densityCsqByPhenotype}}{
    A \code{ggplot} object representing the density of a given VEP prediction
    field (coerced to \code{numeric}) in each level of a given phenotype.}
\item{\code{densityCsqInPhenoLevel}}{
    Similarly, a \code{ggplot} object, for a single level of a given
    phenotype.}
}
}
\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\linkS4class{VCF}},
\code{\link{parseCSQToGRanges,VCF-method}},
\code{\linkS4class{GRanges}},
\code{\linkS4class{BiocParallelParam}},
\code{\link{ensemblVEP}},
and \code{\linkS4class{DataFrame}}.
}

\examples{
# Example data ----

# VCF file
extdata <- file.path(system.file(package = "tSVE"), "extdata")
vcfFile <- file.path(extdata, "moderate.vcf")

# Phenotype file
phenoFile <- file.path(extdata, "moderate_pheno.txt")

# Parameters
tparam <- tSVEParam(
    genos = list(
        REF = "0|0",
        HET = c("0|1", "1|0"),
        ALT = "1|1"))

vcf <- preprocessVariants(
    file = vcfFile, param = tparam, phenos = phenoFile)

# Example usage ----

densityCsqInPhenoLevel(
    level = "EUR", vcf = vcf, phenoCol = "super_pop", csqCol = "CADD_PHRED",
    param = tparam)

densityCsqByPhenotype(
    vcf = vcf, phenoCol = "super_pop", csqCol = "CADD_PHRED",
    param = tparam, facet = "Feature", plot = "TRUE") +
    ggplot2::theme(
        legend.text = ggplot2::element_text(size = ggplot2::rel(.5))
    )
}

\keyword{methods}

