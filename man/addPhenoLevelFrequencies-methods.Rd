\name{addPhenoLevelFrequencies-methods}

\alias{addPhenoLevelFrequencies}
\alias{addPhenoLevelFrequencies,ExpandedVCF-method}

\docType{methods}
\alias{addPhenoLevelFrequencies-methods}
\title{
Genotypes and allele frequencies for a given phenotype level
}

\description{

Adds genotypes counts
(reference homozygote, heterozygote, and alternate homozygote)
and allele frequencies (alternate and minor) calculated
in a group of samples associated with a given level of a given phenotype
as INFO fields in an \code{ExpandedVCF} object.
Multiple genotypes can be counted toward a single
frequency (\emph{e.g.} combined \code{c("0/0", "0|0")} for
homozygote reference genotypes).
}

\usage{
\S4method{addPhenoLevelFrequencies}{ExpandedVCF}(
    vcf, pheno, level, force = FALSE)
}

\arguments{
    \item{vcf}{
\code{ExpandedVCF} object.

\code{metadata(vcf)[["TVTBparam"]]} must contain a
\code{\linkS4class{TVTBparam}} object.
}
    \item{pheno}{
Phenotype in \code{colnames(colData(vcf))}.
}
    \item{level}{
Phenotype level in \code{colData(vcf)[,pheno]}.
}
    \item{force}{
If \code{TRUE}, INFO fields header and data are overwritten.

If \code{FALSE}, an error is thrown if any field already exists.
}
}

\details{
Genotypes and frequencies are calculated within the groups of samples
associated with the given level of the given phenotype,
and stored in INFO fields
named according to settings stored in
\code{metadata(vcf)[["TVTBparam"]]} (see below).

The \code{\linkS4class{TVTBparam}} object controls
the key suffix of INFO fields as follows:
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
\code{\link{addOverallFrequencies,ExpandedVCF-method}},
\code{\linkS4class{VCF}},
and \code{\linkS4class{TVTBparam}}.
}

\examples{
# Example data ----

# VCF file
vcfFile <- system.file("extdata", "moderate.vcf", package = "TVTB")

# Phenotype file
phenoFile <- system.file("extdata", "moderate_pheno.txt", package = "TVTB")
phenotypes <- S4Vectors::DataFrame(read.table(phenoFile, TRUE, row.names = 1))

# TVTB parameters
tparam <- TVTBparam(Genotypes("0|0", c("0|1", "1|0"), "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(
    vcfFile, param = tparam, colData = phenotypes)
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)


# Example usage ----

vcf <- addPhenoLevelFrequencies(vcf, "super_pop", "AFR")
}
