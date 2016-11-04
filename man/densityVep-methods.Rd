\name{densityVep-methods}

\alias{densityVep}

\alias{densityVepByPhenotype}
\alias{densityVepInPhenoLevel}

\alias{densityVepByPhenotype,ExpandedVCF,TVTBparam-method}
\alias{densityVepByPhenotype,ExpandedVCF,missing-method}

\alias{densityVepInPhenoLevel,ExpandedVCF,TVTBparam-method}
\alias{densityVepInPhenoLevel,ExpandedVCF,missing-method}

\docType{methods}
\alias{densityVep-methods}
\alias{densityVepByPhenotype-methods}
\alias{densityVepInPhenoLevel-methods}
\title{Density of VEP predictions in phenotypes}
\description{
\code{densityVepByPhenotype} plots the density of a given VEP prediction
    field (coerced to \code{numeric}) for variants observed (uniquely) in each
    level of a given phenotype.

\code{densityVepInPhenoLevel} considers only variants observed in samples
    associated with a given phenotype level, and plots the corresponding
    density for a given VEP prediction.
}

\usage{
\S4method{densityVepByPhenotype}{ExpandedVCF,TVTBparam}(
        vcf, phenoCol, vepCol, param, ..., filter = VcfFilterRules(),
        unique = FALSE, facet = NULL, plot = FALSE, popFreq = FALSE)

\S4method{densityVepByPhenotype}{ExpandedVCF,missing}(
        vcf, phenoCol, vepCol, alts, param = NULL, ...,
        filter = VcfFilterRules(),
        unique = FALSE, facet = NULL, plot = FALSE, popFreq = FALSE)

\S4method{densityVepInPhenoLevel}{ExpandedVCF,TVTBparam}(
        level, vcf, phenoCol, vepCol, param, ...,
        filter = VcfFilterRules(),
        unique = FALSE, facet = NULL, plot = FALSE, popFreq = FALSE)

\S4method{densityVepInPhenoLevel}{ExpandedVCF,missing}(
        level, vcf, phenoCol, vepCol, alts, param = NULL, ...,
        filter = VcfFilterRules(),
        unique = FALSE, facet = NULL, plot = FALSE, popFreq = FALSE)
}

\arguments{
    \item{vcf}{
\code{ExpandedVCF} object.
}
    \item{phenoCol}{
Name of a column in \code{pheno}.
}
    \item{vepCol}{
Name of a metadata column in \code{ensemblVEP::parseCSQToGRanges(vcf)}.
}
    \item{param}{
\code{TVTBparam} object that contains recurrent parameters.
}
    \item{alts}{
Character vector of alternate genotypes (both heterozygote and homozygote).
}
    \item{level}{
Phenotype level; only variants observed in at least one sample will be
considered.
}
    \item{...}{
Additional arguments used to override \code{TVTBparam} slots.
}
    \item{filter}{
Filter rules of class \code{VcfFixedRules}, \code{VcfInfoRules},
\code{VcfVepRules}, or \code{VcfFilterRules}.
The \code{S4Vectors::subsetByFilter} applies active filter rules
to the \code{vcf} argument prior to summarisation.
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
Should be set to \code{TRUE} if \code{vepCol} is a VEP key representing a
population frequency of the form \code{allele:frequency}.
If set to \code{TRUE}, only \code{frequency} will be extracted and typed as
\code{numeric}.
}
}

\value{
\describe{
\item{\code{densityVepByPhenotype}}{
    A \code{ggplot} object representing the density of a given VEP prediction
    field (coerced to \code{numeric}) in each level of a given phenotype.}
\item{\code{densityVepInPhenoLevel}}{
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

# Pre-process variants
vcf <- VariantAnnotation::readVcf(file = vcfFile)
colData(vcf) <- phenotypes
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)

# Example usage ----

densityVepInPhenoLevel(
    level = "EUR", vcf = vcf, phenoCol = "super_pop",
    vepCol = "CADD_PHRED", param = tparam)

densityVepByPhenotype(
    vcf = vcf, phenoCol = "super_pop", vepCol = "CADD_PHRED",
    param = tparam, facet = "Feature", plot = "TRUE") +
    ggplot2::theme(
        legend.text = ggplot2::element_text(size = ggplot2::rel(.5))
    )
}

\keyword{methods}

