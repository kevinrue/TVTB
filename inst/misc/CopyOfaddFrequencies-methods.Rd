\name{addFrequencies-methods}

\alias{addFrequencies}
\alias{addFrequencies,ExpandedVCF,missing,tSVEParam-method}
\alias{addFrequencies,ExpandedVCF,missing,missing-method}

\alias{addFrequencies,ExpandedVCF,character,tSVEParam-method}
\alias{addFrequencies,ExpandedVCF,character,missing-method}

\alias{addFrequencies,ExpandedVCF,list,tSVEParam-method}
\alias{addFrequencies,ExpandedVCF,list,missing-method}

\docType{methods}
\alias{addFrequencies-methods}
\title{
Group-level genotypes counts and allele frequencies
}
\description{

\code{addFrequencies()} adds, for each level of a phenotype,
    group-level genotypes counts
    (reference homozygote, heterozygote, and alternate homozygote) using
    \code{countGenos()}, and allele frequencies (alternate and minor)
    as \emph{INFO} fields in an \code{ExpandedVCF} object. The function
    requires vectors of valid genotypes expected in the genotype matrix.

\code{addOverallFrequencies()} adds overall genotypes counts
    (reference homozygote, heterozygote, and alternate homozygote) using
    \code{countGenos()}, and overall allele frequencies (alternate and minor)
    as \emph{INFO} fields in an \code{ExpandedVCF} object. The function
    requires vectors of valid genotypes expected in the genotype matrix.

\code{addPhenoLevelFrequencies()} adds, for a given level of a given
    phenotype, group-level genotypes counts and allele frequencies
    as \emph{INFO} fields in a n\code{ExpandedVCF} object. The function
    requires vectors of valid genotypes expected in the genotype matrix.
}
\usage{
addOverallFrequencies(
    vcf, ref, het, alt)

addPhenotypeFrequencies(
    vcf, field, ref, het, alt)

addPhenotypeLevelFrequencies(
    vcf, field, level, ref, het, alt)
}

\arguments{
    \item{vcf}{
\code{ExpandedVCF} object.
}
    \item{field}{
Name of a column in \code{pheno}.
}
    \item{level}{
Phenotype level present in \code{pheno[,field]}.
}
    \item{ref}{
Character vector of reference homozygote genotypes (\emph{e.g.}
    \code{c("0|0", "0/0")}).
}
    \item{het}{
Character vector of heterozygote genotypes (\emph{e.g.}
    \code{c("0|1", "0/1")}).
}
    \item{alt}{
Character vector of alternate homozygote genotypes (\emph{e.g.}
    \code{c("1|1", "1/1")}).
}
}

\details{
\code{addOverallFrequencies()} adds five INFO fields to the VCF object:
\describe{
    \item{\code{REF}}{Count of reference homozygote genotypes.}
    \item{\code{HET}}{Count of heterozygote genotypes.}
    \item{\code{ALT}}{Count of alternate homozygote genotypes.}
    \item{\code{AAF}}{Alternate allele frequency.}
    \item{\code{MAF}}{Minor allele frequency}
}

}

\value{
\describe{
\item{add*Frequencies}{\code{ExpandedVCF} object including additional
    \emph{INFO} fields for genotype counts and allele frequencies.
    See note.}
}
}

\author{
Kevin Rue-Albrecht
}

\note{
Currently, the \code{add*Frequencies} functions do not update existing
    \emph{INFO} fields.
    If any \emph{INFO} field already exists, the functions will return the
    original \code{ExpandedVCF} object.
    In the future, a \code{force} argument might be introduced to allow the
    functions to overwrite existing \emph{INFO} fields, as well as their
    corresponding \code{header} field.
}

\seealso{
\code{\link{addOverallFrequencies,ExpandedVCF,tSVEParam-method}},
\code{\link{addPhenoLevelFrequencies,ExpandedVCF,tSVEParam-method}},
\code{\linkS4class{VCF}}
and \code{\linkS4class{BiocParallelParam}}.
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

# Example usage ----

addOverallFrequencies(
    vcf = eVcf, ref = c("0|0"), het = c("0|1"), alt = c("1|1"))

addPhenotypeFrequencies(
    vcf = eVcf, field = "super_pop",
    ref = c("0|0"), het = c("0|1"), alt = c("1|1"))

addPhenotypeLevelFrequencies(
    vcf = eVcf, level = "AFR", field = "super_pop",
    ref = c("0|0"), het = c("0|1"), alt = c("1|1"))
}
