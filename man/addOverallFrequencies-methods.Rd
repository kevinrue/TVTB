\name{addOverallFrequencies-methods}

\alias{addOverallFrequencies}
\alias{addOverallFrequencies,ExpandedVCF,tSVEParam-method}
\alias{addOverallFrequencies,ExpandedVCF,missing-method}

\docType{methods}
\alias{addOverallFrequencies-methods}
\title{
Overall genotypes counts and allele frequencies
}

\description{

Adds dataset-wide genotypes counts
(reference homozygote, heterozygote, and alternate homozygote)
and allele frequencies (alternate and minor)
as \emph{INFO} fields in an \code{ExpandedVCF} object.
Multiple genotypes can be counted toward a single
frequency (\emph{e.g.} combined \code{c("0/0", "0|0")} for
homozygote reference genotypes).
}

\usage{
\S4method{addOverallFrequencies}{ExpandedVCF,tSVEParam}(
    vcf, param, ..., force = FALSE)

\S4method{addOverallFrequencies}{ExpandedVCF,missing}(
    vcf, ref, het, alt, ..., force = FALSE)
}

\arguments{
    \item{vcf}{
\code{ExpandedVCF} object.
}
    \item{param}{
\code{tSVEParam} object that contains recurrent parameters.
}
    \item{ref}{
\code{character} vector of valid reference homozygote genotypes (\emph{e.g.}
    \code{c("0|0", "0/0")}).
}
    \item{het}{
\code{character} vector of valid heterozygote genotypes (\emph{e.g.}
    \code{c("0|1", "0/1")}).
}
    \item{alt}{
\code{character} vector of valid alternate homozygote genotypes (\emph{e.g.}
    \code{c("1|1", "1/1")}).
}
    \item{...}{
Additional arguments used to override \code{tSVEParam} slots.
}
    \item{force}{
If \code{TRUE}, \emph{INFO} fields header and data are overwritten.

If \code{FALSE}, an error is thrown if any field already exists.
}
}

\details{
Genotypes and frequencies are calculated across all the
samples in the \code{ExpandedVCF} object, and stored in \emph{INFO} fields
named according to settings stored in the \code{tSVEParam} object (see below).

The \code{param} object controls the name of \emph{INFO} fields as
follows:
\describe{
    \item{\code{names(hRef(param))}}{
        Count of reference homozygote genotypes.}
    \item{\code{names(het(param))}}{
        Count of heterozygote genotypes.}
    \item{\code{names(hAlt(param))}}{
        Count of alternate homozygote genotypes.}
    \item{\code{aaf(param)}}{
        Alternate allele frequency.}
    \item{\code{maf(param)}}{
        Minor allele frequency}
}
Alternatively, if \code{param} is not provided, the following default
values are used respectively:
\code{REF}, \code{HET}, \code{ALT}, \code{AAF}, \code{MAF}.

}

\value{
\code{ExpandedVCF} object including additional
\emph{INFO} fields for genotype counts and allele frequencies.
See \emph{Details}.
}

\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\link{addFrequencies,ExpandedVCF,list,tSVEParam-method}},
\code{\link{addPhenoLevelFrequencies,ExpandedVCF,tSVEParam-method}},
and \code{\linkS4class{VCF}}.
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

svp <- VariantAnnotation::ScanVcfParam(
    geno = "GT",
    info = "CSQ",
    samples = rownames(phenotypes),
    which = bedRegions)
tf <- Rsamtools::TabixFile(vcfFile)
vcf <- VariantAnnotation::readVcf(file = tf, param = svp)
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

addOverallFrequencies(
    vcf = eVcf,
    param = tparam)
}
