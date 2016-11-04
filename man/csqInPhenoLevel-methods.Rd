\name{csqInPhenoLevel-methods}

\alias{csqInPhenoLevel}
\alias{csqInPhenoLevel,ExpandedVCF,tSVEParam-method}
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
\S4method{csqInPhenoLevel}{ExpandedVCF,tSVEParam}(
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
\code{tSVEParam} object that contains recurrent parameters.
}
    \item{alts}{
\code{character} vector of alternate genotypes
(both heterozygote and homozygote).
}
    \item{...}{
Additional arguments used to override \code{tSVEParam} slots.
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

csqInPhenoLevel(
    vcf = eVcf, phenoCol = "super_pop", level = "AFR",
    csqCol = "CADD_PHRED", param = tparam)
}

\keyword{methods}

