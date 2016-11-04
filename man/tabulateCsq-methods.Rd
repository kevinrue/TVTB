\name{tabulateCsq-methods}

\alias{tabulateCsq}

\alias{tabulateCsqByPhenotype}
\alias{tabulateCsqInPhenoLevel}

\alias{tabulateCsqByPhenotype,ExpandedVCF,tSVEParam-method}
\alias{tabulateCsqByPhenotype,ExpandedVCF,missing-method}

\alias{tabulateCsqInPhenoLevel,ExpandedVCF,tSVEParam-method}
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
\S4method{tabulateCsqByPhenotype}{ExpandedVCF,tSVEParam}(
        vcf, phenoCol, csqCol, param, ...,
        unique = FALSE, facet = NULL, plot = FALSE, percentage = FALSE)

\S4method{tabulateCsqByPhenotype}{ExpandedVCF,missing}(
        vcf, phenoCol, csqCol, alts, param = NULL, ...,
        unique = FALSE, facet = NULL, plot = FALSE, percentage = FALSE)

## tabulateCsqInPhenoLevel
\S4method{tabulateCsqInPhenoLevel}{ExpandedVCF,tSVEParam}(
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
library(VariantAnnotation)

# Example data ----

extdata <- file.path(system.file(package = "tSVE"), "extdata")
vcfFile <- file.path(extdata, "chr15.phase3_integrated.vcf.gz")
phenoFile <- file.path(extdata, "integrated_samples.txt")
bedRegions <- GenomicRanges::GRanges(
    seqnames = "15",
    ranges = IRanges::IRanges(start = 48420E3, end = 48421E3))

# Import phenotype data ----

phenotypes <- S4Vectors::DataFrame(read.table(
    file = phenoFile,
    header = TRUE,
    row.names = 1))

# Import variants ----

svp <- ScanVcfParam(
    geno = "GT",
    info = "CSQ",
    samples = rownames(phenotypes),
    which = bedRegions)
tf <- Rsamtools::TabixFile(vcfFile)
vcf <- readVcf(file = tf, param = svp)
# Separate multi-allelic records into bi-allelic records
eVcf <- expand(x = vcf, row.names = TRUE)
# Disambiguate row.names from multi-allelic records
rownames(eVcf) <- paste(rownames(eVcf), mcols(eVcf)[,"ALT"], sep = "_")
# Add some phenotypes information necessary for the demo
colData(eVcf) <- phenotypes

# Define genotypes ----
tparam <- tSVEParam(genos = list(
    REF = c("0|0"),
    HET = c("0|1","1|0"),
    ALT = c("1|1")))

# Example usage ----

# Count of each "Consequence" in "super_pop"="AFR"
tabulateCsqInPhenoLevel(
    level = "AFR", vcf = eVcf, phenoCol = "super_pop",
    csqCol = "Consequence", param = tparam,
    facet = "Feature")

# Count of each "Consequence" in each "super_pop" (barplot)
tabulateCsqByPhenotype(
    vcf = eVcf, phenoCol = "super_pop",
    csqCol = "Consequence", param = tparam,
    facet = "Feature", plot = TRUE) +
    ggplot2::theme(
        legend.text = ggplot2::element_text(size = ggplot2::rel(.5))
    )

# Count of each "Consequence" in each "super_pop" (data.frame)
as.data.frame(table(
    tabulateCsqByPhenotype(
        vcf = eVcf, phenoCol = "super_pop",
        csqCol = "Consequence", param = tparam,
        facet = "Feature", plot = FALSE)
))
}

\keyword{methods}

