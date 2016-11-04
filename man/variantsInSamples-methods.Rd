\name{variantsInSamples-methods}

\alias{variantsInSamples}
\alias{variantsInSamples,ExpandedVCF,TVTBparam-method}

\docType{methods}
\alias{variantsInSamples-methods}
\title{Identify variants observed in samples}
\description{
Identifies variants observed (uniquely) in at least one sample of a given
group.
}

\usage{
\S4method{variantsInSamples}{ExpandedVCF,TVTBparam}(
    vcf, param, samples = 1:ncol(vcf), ..., unique = FALSE)

\S4method{variantsInSamples}{ExpandedVCF,TVTBparam}(
    vcf, param, samples = 1:ncol(vcf), ..., unique = FALSE)
}

\arguments{
    \item{vcf}{
\code{ExpandedVCF} object.
}
    \item{samples}{
\code{integer}, \code{numeric} or \code{character} vector
indicating samples to consider in \code{VariantAnnotation::geno(vcf)}.
If not specified, all samples are considered.
}
    \item{param}{
\code{TVTBparam} object that contains recurrent parameters.
}
    \item{...}{
Additional arguments used to override \code{TVTBparam} slots (\emph{e.g.}
\code{het}, \code{alt}).
}
    \item{unique}{
If \code{TRUE}, consider only variants \emph{unique} to the phenotype level
(\emph{i.e.} not seen in any other phenotype level).
}
}

\value{
Am \code{integer} vector of indices indicating which variants are observed
in at least one non-reference genotype in the given group of samples.
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

variantsInSamples(
    vcf = vcf,
    samples = which(SummarizedExperiment::colData(vcf)[,"super_pop"]=="EUR"),
    param = tparam)
}

\keyword{methods}

