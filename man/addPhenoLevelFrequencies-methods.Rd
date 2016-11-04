\name{addPhenoLevelFrequencies-methods}

\alias{addPhenoLevelFrequencies}
\alias{addPhenoLevelFrequencies,ExpandedVCF,tSVEParam-method}
\alias{addPhenoLevelFrequencies,ExpandedVCF,missing-method}

\docType{methods}
\alias{addPhenoLevelFrequencies-methods}
\title{
Genotypes and allele frequencies for a given phenotype level
}

\description{

Adds genotypes counts
(reference homozygote, heterozygote, and alternate homozygote)
and allele frequencies (alternate and minor)
in a group of samples associated with a given level of a given phenotype
as \emph{INFO} fields in an \code{ExpandedVCF} object.
Multiple genotypes can be counted toward a single
frequency (\emph{e.g.} combined \code{c("0/0", "0|0")} for
homozygote reference genotypes).
}

\usage{
\S4method{addPhenoLevelFrequencies}{ExpandedVCF,tSVEParam}(
    vcf, pheno, level, param, ..., force = FALSE)

\S4method{addPhenoLevelFrequencies}{ExpandedVCF,missing}(
    vcf, pheno, level, ref, het, alt, ..., force = FALSE)
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
Genotypes and frequencies are calculated within the groups of samples
associated with the given level of the given phenotype,
and stored in \emph{INFO} fields
named according to settings stored in the \code{tSVEParam} object (see below).

The \code{param} object controls the suffix of \emph{INFO} fields as
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
\code{\link{addOverallFrequencies,ExpandedVCF,tSVEParam-method}},
and \code{\linkS4class{VCF}}.
}

\examples{
# Example data ----

# VCF file
extdata <- file.path(system.file(package = "tSVE"), "extdata")
vcfFile <- file.path(extdata, "moderate.vcf")

# Phenotype file
phenoFile <- file.path(extdata, "moderate_pheno.txt")

# Parameters
tparam <- tSVEParam(
    genos = list(
        REF = "0|0",
        HET = c("0|1", "1|0"),
        ALT = "1|1"))

vcf <- preprocessVariants(
    file = vcfFile, param = tparam, phenos = phenoFile)

# Example usage ----

addPhenoLevelFrequencies(
    vcf = vcf,
    pheno = "super_pop",
    level = "AFR",
    param = tparam)
}
