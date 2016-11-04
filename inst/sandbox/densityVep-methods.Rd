\name{densityVep-methods}

\alias{densityVep}

\alias{densityVepByPhenotype}
\alias{densityVepInPhenoLevel}

\alias{densityVepByPhenotype,ExpandedVCF,TVTBparam-method}
\alias{densityVepInPhenoLevel,ExpandedVCF,TVTBparam-method}

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
        unique = FALSE, facet = NULL, plot = FALSE, pattern = NULL,
        layer = "density+dotplot")

\S4method{densityVepInPhenoLevel}{ExpandedVCF,TVTBparam}(
        level, vcf, phenoCol, vepCol, param, ...,
        filter = VcfFilterRules(),
        unique = FALSE, facet = NULL, plot = FALSE, pattern = NULL,
        layer = "density+dotplot")
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
\code{TVTBparam} object that contains recurrent parameters (\emph{e.g.}
\code{het}, \code{alt}).
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
as teh Ensembl VEP typically produces one prediction per variant per feature.
}
    \item{plot}{
If \code{TRUE}, return a \code{ggplot} object.
}
    \item{pattern}{
If \code{NULL}, \code{vepCol} is assumed to be immediately coercible
to \code{numeric}.
Otherwise, a regex pattern (passed to \code{gsub})
can be provided to extract a substring that must
be coercible to \code{numeric}.

For instance, the patterns
\code{"[[:alpha:]]*:(.*)"} and \code{".*:(.*)"} are both valid
patterns to extract the population allele frequency \code{0.008}
from a value stored in the form \code{"A:0.008"}
(more generally, \code{allele:frequency}).
A single pair of round brackets \code{()} is expected
in the pattern,
to capture the substring that will be coerced to \code{numeric}.
}
    \item{layer}{
A \code{character} value of \code{"+"}-separated \code{geom} types to apply as
layers if \code{plot=TRUE}. Currently, only \code{density} and \code{dotplot}
are accepted (\emph{e.g.} \code{layer="density+dotplot"}).
Note that an empty string can be used to obtain the \code{ggplot} object
without any layer.
}
}

\value{
\describe{
\item{\code{densityVepByPhenotype}}{
    \code{plot=TRUE}: A \code{ggplot} object representing the density of a
    given VEP prediction field (coerced to \code{numeric} after optional
    pattern extraction) in each level of a given phenotype.

    \code{plot=FALSE}: A \code{data.frame} in long format including one row per
    variant seen in at least one sample (heterozygote or alternate homozygote)
    associated with each phenotype level.
    The \code{data.frame} contains at least one column for the VEP prediction
    value coerced to \code{numeric} after the optional pattern extraction,
    and one column for the phenotype level. An additional column
    containing another VEP prediction field may be added using the \code{facet}
    argument.}

\item{\code{densityVepInPhenoLevel}}{
    Similarly to \code{tabulateVepByPhenotype} above, either a \code{ggplot}
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
vcfFile <- system.file("extdata", "moderate.vcf", package = "TVTB")

# Phenotype file
phenoFile <- system.file("extdata", "moderate_pheno.txt", package = "TVTB")
phenotypes <- S4Vectors::DataFrame(
    read.table(file = phenoFile, header = TRUE, row.names = 1))

# TVTB parameters
tparam <- TVTBparam(Genotypes("0|0", c("0|1", "1|0"), "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(vcfFile)
colData(vcf) <- phenotypes
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)


# Example usage ----

densityVepInPhenoLevel("EUR", vcf, "super_pop", "CADD_PHRED", tparam)

densityVepByPhenotype(
    vcf, "super_pop", "CADD_PHRED", tparam, facet = "Feature", plot = "TRUE"
    ) +
    ggplot2::theme(
        legend.text = ggplot2::element_text(size = ggplot2::rel(.5))
    )
}

\keyword{methods}
