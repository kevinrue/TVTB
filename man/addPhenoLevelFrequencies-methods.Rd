\name{addPhenoLevelFrequencies-methods}

\alias{addPhenoLevelFrequencies}
\alias{addPhenoLevelFrequencies,ExpandedVCF,TVTBparam-method}

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
\S4method{addPhenoLevelFrequencies}{ExpandedVCF,TVTBparam}(
    vcf, pheno, level, param, ..., force = FALSE)
}

\arguments{
    \item{vcf}{
\code{ExpandedVCF} object.
}
    \item{pheno}{
Phenotype in \code{colnames(colData(vcf))}.
}
    \item{level}{
Phenotype level in \code{colData(vcf)[,pheno]}.
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
Genotypes and frequencies are calculated within the groups of samples
associated with the given level of the given phenotype,
and stored in INFO fields
named according to settings stored in the \code{TVTBparam} object (see below).

The \code{param} object controls the key suffix of INFO fields as
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
\code{\link{addOverallFrequencies,ExpandedVCF,TVTBparam-method}},
and \code{\linkS4class{VCF}}.
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

# Pre-process variants
vcf <- VariantAnnotation::readVcf(file = vcfFile)
colData(vcf) <- phenotypes
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)

# Example usage ----

addPhenoLevelFrequencies(
    vcf = vcf,
    pheno = "super_pop",
    level = "AFR",
    param = tparam)
}
