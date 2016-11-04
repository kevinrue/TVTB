\name{variantsInSamples-methods}

\alias{variantsInSamples}
\alias{variantsInSamples,ExpandedVCF,numeric,tSVEParam-method}
\alias{variantsInSamples,ExpandedVCF,numeric,missing-method}
\alias{variantsInSamples,ExpandedVCF,character,tSVEParam-method}
\alias{variantsInSamples,ExpandedVCF,character,missing-method}

\docType{methods}
\alias{variantsInSamples-methods}
\title{Identify variants observed in samples}
\description{
Identifies variants observed (uniquely) in at least one sample of a given
group.
}

\usage{
\S4method{variantsInSamples}{ExpandedVCF,numeric,tSVEParam}(
    vcf, samples, param, ..., unique = FALSE)

\S4method{variantsInSamples}{ExpandedVCF,character,tSVEParam}(
    vcf, samples, param, ..., unique = FALSE)

\S4method{variantsInSamples}{ExpandedVCF,numeric,missing}(
    vcf, samples, alts, unique = FALSE)

\S4method{variantsInSamples}{ExpandedVCF,character,missing}(
    vcf, samples, alts, unique = FALSE)
}

\arguments{
    \item{vcf}{
\code{ExpandedVCF} object.
}
    \item{samples}{
\code{numeric} or \code{character} vector of indicating samples to consider
in \code{VariantAnnotation::geno(vcf)}.
}
    \item{param}{
\code{tSVEParam} object that contains recurrent parameters.
}
    \item{alts}{
\code{character} vector of alternate genotypes (both heterozygote and
homozygote).
}
    \item{...}{
Additional arguments used to override \code{tSVEParam} slots.
}
    \item{unique}{
If \code{TRUE}, consider only variants \emph{unique} to the phenotype level
(\emph{i.e.} not seen in any other phenotype level).
}
}

\value{
A numeric vector of indices indicating which variants are observed in at least
one non-reference genotype in the given group.
}
\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\linkS4class{VCF}},
\code{\link{geno,VCF-method}}
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

# Define genotypes
tparam <- tSVEParam(genos = list(
    REF = c("0|0"),
    HET = c("0|1","1|0"),
    ALT = c("1|1")))

# Example usage ----

variantsInSamples(
    vcf = eVcf,
    samples = which(colData(eVcf)[,"super_pop"] == "EUR"),
    param = tparam)
}

\keyword{methods}

