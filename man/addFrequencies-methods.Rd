\name{addFrequencies-methods}

\alias{addFrequencies}

\alias{addFrequencies,ExpandedVCF,list,TVTBparam-method}
\alias{addFrequencies,ExpandedVCF,character,TVTBparam-method}
\alias{addFrequencies,ExpandedVCF,missing,TVTBparam-method}

\docType{methods}
\alias{addFrequencies-methods}
\title{
Group-level genotypes counts and allele frequencies
}

\description{
Adds genotypes counts
(reference homozygote, heterozygote, and alternate homozygote)
and allele frequencies (alternate and minor) 
as INFO fields in an \code{ExpandedVCF} object.
Counts and frequencies may be calculated overall
(\emph{i.e.} across all samples), or within groups of samples
(\emph{i.e.} within phenotype levels).
Multiple genotypes can be counted toward a single
frequency (\emph{e.g.} combined \code{c("0/0", "0|0")} for
homozygote reference genotypes).
}

\usage{
\S4method{addFrequencies}{ExpandedVCF,list,TVTBparam}(
        vcf, phenos, param, ..., force = FALSE)

\S4method{addFrequencies}{ExpandedVCF,character,TVTBparam}(
        vcf, phenos, param, ..., force = FALSE)

\S4method{addFrequencies}{ExpandedVCF,missing,TVTBparam}(
        vcf, param, ..., force = FALSE)
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
See \emph{Details} below.
}
    \item{param}{
\code{TVTBparam} object that contains recurrent parameters.
}
    \item{...}{
Additional arguments used to override \code{TVTBparam} slots (\emph{e.g.}
\code{ref}, \code{het}, \code{alt}).
}
    \item{force}{
If \code{TRUE}, INFO fields header and data are overwritten
with a message, if present.

If \code{FALSE}, an error is thrown if any field already exists.
}
}

\details{
The \code{phenos} argument is central to control the behaviour of this method.

If \code{phenos=NULL}, genotypes and frequencies are calculated across all the
samples in the \code{ExpandedVCF} object, and stored in INFO fields
named according to settings stored in the \code{TVTBparam} object (see below).

If \code{phenos} is a \code{character} vector of phenotypes present in
\code{colnames(colData(vcf))}, counts and frequencies are calculated for each
level of those phenotypes, and stored in INFO fields prefixed with
"<phenotype>_<level>_" and suffixed with the settings stored in the
\code{param} object (see below).

Finally, if \code{phenos} is a named \code{list}, \code{names} must be
phenotypes present in \code{colnames(colData(vcf))}, and values must be levels
of those phenotypes. In this case, counts and frequencies are calculated for
the given levels of the given phenotypes, and stored in INFO fields as
described above.

The \code{param} object controls the key (suffix) of INFO fields as
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
\code{\link{addOverallFrequencies,ExpandedVCF,TVTBparam-method}},
\code{\link{addPhenoLevelFrequencies,ExpandedVCF,TVTBparam-method}},
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

addFrequencies(
    vcf = vcf,
    phenos = list(super_pop = "AFR"),
    param = tparam)
}
