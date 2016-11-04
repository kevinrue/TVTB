\name{addOverallFrequencies-methods}

\alias{addOverallFrequencies}
\alias{addOverallFrequencies,ExpandedVCF-method}

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
\S4method{addOverallFrequencies}{ExpandedVCF}(vcf, force = FALSE)
}

\arguments{
    \item{vcf}{
\code{ExpandedVCF} object.

\code{metadata(vcf)[["TVTBparam"]]} must contain a
\code{\linkS4class{TVTBparam}} object.
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
    \item{\code{names(ref(param))}}{
        Count of reference homozygote genotypes.}
    \item{\code{names(het(param))}}{
        Count of heterozygote genotypes.}
    \item{\code{names(alt(param))}}{
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

\section{Warning}{
A \code{warning} message is issued if genotypes are not fully defined in the
\code{TVTBparam}.
}

\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\link{addFrequencies,ExpandedVCF,list-method}},
\code{\link{addPhenoLevelFrequencies,ExpandedVCF-method}},
and \code{\linkS4class{VCF}}.
}

\examples{
# Example data ----

# VCF file
vcfFile <- system.file("extdata", "moderate.vcf", package = "TVTB")

# TVTB parameters
tparam <- TVTBparam(Genotypes("0|0", c("0|1", "1|0"), "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(vcfFile, param = tparam)
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)


# Example usage ----

vcf <- addOverallFrequencies(vcf, tparam)
}
