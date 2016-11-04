\name{addOverallFrequencies-methods}

\alias{addOverallFrequencies}
\alias{addOverallFrequencies,ExpandedVCF,TVTBparam-method}

\docType{methods}
\alias{addOverallFrequencies-methods}
\title{
Overall genotypes counts and allele frequencies
}

\description{

Adds dataset-wide genotypes counts
(reference homozygote, heterozygote, and alternate homozygote)
and allele frequencies (alternate and minor)
as INFO fields in an \code{ExpandedVCF} object.
Counts and frequencies may be calculated across all samples.
Multiple genotypes can be counted toward a single
frequency (\emph{e.g.} combined \code{c("0/0", "0|0")} for
homozygote reference genotypes).
}

\usage{
\S4method{addOverallFrequencies}{ExpandedVCF,TVTBparam}(
    vcf, param, ..., force = FALSE)
}

\arguments{
    \item{vcf}{
\code{ExpandedVCF} object.
}
    \item{param}{
\code{TVTBparam} object that contains recurrent parameters.
}
    \item{...}{
Additional arguments used to override \code{TVTBparam} slots (\emph{e.g.}
\code{ref}, \code{het}, \code{alt}).
}
    \item{force}{
If \code{TRUE}, INFO fields header and data are overwritten.

If \code{FALSE}, an error is thrown if any field already exists.
}
}

\details{
Genotypes and frequencies are calculated across all the
samples in the \code{ExpandedVCF} object, and stored in INFO fields
named according to settings stored in the \code{TVTBparam} object (see below).

The \code{param} object controls the key of INFO fields as
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
}

\value{
\code{ExpandedVCF} object including additional
INFO fields for genotype counts and allele frequencies.
See \emph{Details}.
}

\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\link{addFrequencies,ExpandedVCF,list,TVTBparam-method}},
\code{\link{addPhenoLevelFrequencies,ExpandedVCF,TVTBparam-method}},
and \code{\linkS4class{VCF}}.
}

\examples{
# Example data ----

# VCF file
extdata <- file.path(system.file(package = "TVTB"), "extdata")
vcfFile <- file.path(extdata, "moderate.vcf")

# TVTB parameters
tparam <- TVTBparam(
    genos = list(
        REF = "0|0",
        HET = c("0|1", "1|0"),
        ALT = "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(file = vcfFile)
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)

# Example usage ----

addOverallFrequencies(
    vcf = vcf,
    param = tparam)
}
