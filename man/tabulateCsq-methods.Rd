\name{tabulateCsq-methods}

\alias{tabulateCsq}

\alias{tabulateCsqByPhenotype}
\alias{tabulateCsqInPhenoLevel}

\alias{tabulateCsqByPhenotype,ExpandedVCF,TVTBparam-method}
\alias{tabulateCsqByPhenotype,ExpandedVCF,missing-method}

\alias{tabulateCsqInPhenoLevel,ExpandedVCF,TVTBparam-method}
\alias{tabulateCsqInPhenoLevel,ExpandedVCF,missing-method}

\docType{methods}
\alias{tabulateCsq-methods}
\alias{tabulateCsqByPhenotype-methods}
\alias{tabulateCsqInPhenoLevel-methods}
\title{Tabulate VEP predictions in phenotypes}
\description{
\code{tabulateCsqByPhenotype} tabulates or plots counts of a given VEP
    prediction field for variants observed (uniquely) in each level of a given
    phenotype.

\code{tabulateCsqInPhenoLevel} considers only variants observed in samples
    associated with a given phenotype level, and tabulates or plots the
    corresponding counts for a given VEP prediction.
}

\usage{
## tabulateCsqByPhenotype
\S4method{tabulateCsqByPhenotype}{ExpandedVCF,TVTBparam}(
        vcf, phenoCol, csqCol, param, ...,
        unique = FALSE, facet = NULL, plot = FALSE, percentage = FALSE)

\S4method{tabulateCsqByPhenotype}{ExpandedVCF,missing}(
        vcf, phenoCol, csqCol, alts, param = NULL, ...,
        unique = FALSE, facet = NULL, plot = FALSE, percentage = FALSE)

## tabulateCsqInPhenoLevel
\S4method{tabulateCsqInPhenoLevel}{ExpandedVCF,TVTBparam}(
        level, vcf, phenoCol, csqCol, param, ...,
        unique = FALSE, facet = NULL, plot = FALSE, percentage = FALSE)

\S4method{tabulateCsqInPhenoLevel}{ExpandedVCF,missing}(
        level, vcf, phenoCol, csqCol, alts, param = NULL, ...,
        unique = FALSE, facet = NULL, plot = FALSE, percentage = FALSE)
}

\arguments{
    \item{vcf}{
\code{ExpandedVCF} object.
}
    \item{phenoCol}{
Name of a column in \code{colData(vcf)}.
}
    \item{csqCol}{
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
    \item{percentage}{
\emph{Only considered if \code{plot=TRUE}}. If \code{TRUE}, the counts of VEP
predictions will be represented as stacked percentage values.
}
}

\value{
\describe{
\item{\code{tabulateCsqByPhenotype}}{
    \code{plot=TRUE}: A \code{ggplot} object representing the counts of
    values for a given VEP prediction field in each level of a given phenotype.

    \code{plot=FALSE}: A \code{data.frame} in long format including one row per
    variant seen in at least one sample (heterozygote or alternate homozygote)
    associated with each phenotype level.
    The \code{data.frame} contains at least one column for the VEP prediction
    value, and one column for the phenotype level. An additional column
    containing another VEP prediction field may be added using the \code{facet}
    argument.}
\item{\code{tabulateCsqInPhenoLevel}}{
    Similarly to \code{tabulateCsqByPhenotype} above, either a \code{ggplot}
    object or a \code{data.frame} in long format, for a single level of
    a given phenotype.}
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

# Parameters
tparam <- TVTBparam(
    genos = list(
        REF = "0|0",
        HET = c("0|1", "1|0"),
        ALT = "1|1"))

vcf <- preprocessVariants(
    file = vcfFile, param = tparam, phenos = phenoFile)

# Example usage ----

# Count of each "Consequence" in "super_pop"="AFR"
tabulateCsqInPhenoLevel(
    level = "AFR", vcf = vcf, phenoCol = "super_pop",
    csqCol = "Consequence", param = tparam,
    facet = "Feature")

# Count of each "Consequence" in each "super_pop" (barplot)
tabulateCsqByPhenotype(
    vcf = vcf, phenoCol = "super_pop",
    csqCol = "Consequence", param = tparam,
    facet = "Feature", plot = TRUE) +
    ggplot2::theme(
        legend.text = ggplot2::element_text(size = ggplot2::rel(.5))
    )

# Count of each "Consequence" in each "super_pop" (data.frame)
as.data.frame(table(
    tabulateCsqByPhenotype(
        vcf = vcf, phenoCol = "super_pop",
        csqCol = "Consequence", param = tparam,
        facet = "Feature", plot = FALSE)
))
}

\keyword{methods}

