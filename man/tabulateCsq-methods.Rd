\name{tabulateVep-methods}

\alias{tabulateVep}

\alias{tabulateVepByPhenotype}
\alias{tabulateVepInPhenoLevel}

\alias{tabulateVepByPhenotype,ExpandedVCF,TVTBparam-method}
\alias{tabulateVepByPhenotype,ExpandedVCF,missing-method}

\alias{tabulateVepInPhenoLevel,ExpandedVCF,TVTBparam-method}
\alias{tabulateVepInPhenoLevel,ExpandedVCF,missing-method}

\docType{methods}
\alias{tabulateVep-methods}
\alias{tabulateVepByPhenotype-methods}
\alias{tabulateVepInPhenoLevel-methods}
\title{Tabulate VEP predictions in phenotypes}
\description{
\code{tabulateVepByPhenotype} tabulates or plots counts of a given VEP
    prediction field for variants observed (uniquely) in each level of a given
    phenotype.

\code{tabulateVepInPhenoLevel} considers only variants observed in samples
    associated with a given phenotype level, and tabulates or plots the
    corresponding counts for a given VEP prediction.
}

\usage{
## tabulateVepByPhenotype
\S4method{tabulateVepByPhenotype}{ExpandedVCF,TVTBparam}(
        vcf, phenoCol, vepCol, param, ...,
        unique = FALSE, facet = NULL, plot = FALSE, percentage = FALSE)

\S4method{tabulateVepByPhenotype}{ExpandedVCF,missing}(
        vcf, phenoCol, vepCol, alts, param = NULL, ...,
        unique = FALSE, facet = NULL, plot = FALSE, percentage = FALSE)

## tabulateVepInPhenoLevel
\S4method{tabulateVepInPhenoLevel}{ExpandedVCF,TVTBparam}(
        level, vcf, phenoCol, vepCol, param, ...,
        unique = FALSE, facet = NULL, plot = FALSE, percentage = FALSE)

\S4method{tabulateVepInPhenoLevel}{ExpandedVCF,missing}(
        level, vcf, phenoCol, vepCol, alts, param = NULL, ...,
        unique = FALSE, facet = NULL, plot = FALSE, percentage = FALSE)
}

\arguments{
    \item{vcf}{
\code{ExpandedVCF} object.
}
    \item{phenoCol}{
Name of a column in \code{colData(vcf)}.
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
\item{\code{tabulateVepByPhenotype}}{
    \code{plot=TRUE}: A \code{ggplot} object representing the counts of
    values for a given VEP prediction field in each level of a given phenotype.

    \code{plot=FALSE}: A \code{data.frame} in long format including one row per
    variant seen in at least one sample (heterozygote or alternate homozygote)
    associated with each phenotype level.
    The \code{data.frame} contains at least one column for the VEP prediction
    value, and one column for the phenotype level. An additional column
    containing another VEP prediction field may be added using the \code{facet}
    argument.}
\item{\code{tabulateVepInPhenoLevel}}{
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

# Count of each "Consequence" in "super_pop"="AFR"
tabulateVepInPhenoLevel(
    level = "AFR", vcf = vcf, phenoCol = "super_pop",
    vepCol = "Consequence", param = tparam,
    facet = "Feature")

# Count of each "Consequence" in each "super_pop" (barplot)
tabulateVepByPhenotype(
    vcf = vcf, phenoCol = "super_pop",
    vepCol = "Consequence", param = tparam,
    facet = "Feature", plot = TRUE) +
    ggplot2::theme(
        legend.text = ggplot2::element_text(size = ggplot2::rel(.5))
    )

# Count of each "Consequence" in each "super_pop" (data.frame)
as.data.frame(table(
    tabulateVepByPhenotype(
        vcf = vcf, phenoCol = "super_pop",
        vepCol = "Consequence", param = tparam,
        facet = "Feature", plot = FALSE)
))
}

\keyword{methods}

