\name{addFrequencies-methods}

\alias{addFrequencies}

\alias{addFrequencies,ExpandedVCF,list,tSVEParam-method}
\alias{addFrequencies,ExpandedVCF,list,missing-method}

\alias{addFrequencies,ExpandedVCF,character,tSVEParam-method}
\alias{addFrequencies,ExpandedVCF,character,missing-method}

\alias{addFrequencies,ExpandedVCF,missing,tSVEParam-method}
\alias{addFrequencies,ExpandedVCF,missing,missing-method}

\docType{methods}
\alias{addFrequencies-methods}
\title{
Group-level genotypes counts and allele frequencies
}

\description{

Adds group-level genotypes counts
(reference homozygote, heterozygote, and alternate homozygote)
and allele frequencies (alternate and minor)
as \emph{INFO} fields in an \code{ExpandedVCF} object.
Multiple genotypes can be counted toward a single
frequency (\emph{e.g.} combined \code{c("0/0", "0|0")} for
homozygote reference genotypes).
}

\usage{
\S4method{addFrequencies}{ExpandedVCF,list,tSVEParam}(
        vcf, phenos, param, ..., force = FALSE)

\S4method{addFrequencies}{ExpandedVCF,list,missing}(
        vcf, ref, het, alt, phenos, ..., force = FALSE)

\S4method{addFrequencies}{ExpandedVCF,character,tSVEParam}(
        vcf, phenos, param, ..., force = FALSE)

\S4method{addFrequencies}{ExpandedVCF,character,missing}(
        vcf, ref, het, alt, phenos, ..., force = FALSE)

\S4method{addFrequencies}{ExpandedVCF,missing,tSVEParam}(
        vcf, param, ..., force = FALSE)

\S4method{addFrequencies}{ExpandedVCF,missing,missing}(
        vcf, ref, het, alt, ..., force = FALSE)
}

\arguments{
    \item{vcf}{
\code{ExpandedVCF} object.
}
    \item{phenos}{
If \code{NULL}, counts and frequencies are calculated across all samples.

Otherwise,
either a \code{character} vector of phenotypes in
\code{colnames(colData(vcf))},
or a named \code{list} in which \code{names} are phenotypes in
\code{colnames(colData(vcf))} and values are \code{character} vectors of
phenotype levels in \code{colData(vcf)[,phenotype]}.
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
The \code{phenos} argument is central to control the behaviour of this method.

If \code{phenos=NULL}, genotypes and frequencies are calculated across all the
samples in the \code{ExpandedVCF} object, and stored in \emph{INFO} fields
named according to settings stored in the \code{tSVEParam} object (see below).

If \code{phenos} is a \code{character} vector of phenotypes present in
\code{colnames(colData(vcf))}, counts and frequencies are calculated for each
level of those phenotypes, and stored in \emph{INFO} fields prefixed with
"<phenotype>_<level>_" and suffixed with the settings stored in the
\code{param} object (see below).

Finally, if \code{phenos} is a named \code{list}, \code{names} must be
phenotypes present in \code{colnames(colData(vcf))}, and values must be levels
of those phenotypes. In this case, counts and frequencies are calculated for
the given levels of the given phenotypes, and stored in \emph{INFO} fields as
described above.

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
\code{\link{addOverallFrequencies,ExpandedVCF,tSVEParam-method}},
\code{\link{addPhenoLevelFrequencies,ExpandedVCF,tSVEParam-method}},
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

addFrequencies(
    vcf = vcf,
    phenos = list(super_pop = "EUR"),
    param = tparam)
}
