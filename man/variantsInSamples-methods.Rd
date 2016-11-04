\name{variantsInSamples-methods}

\alias{variantsInSamples}
\alias{variantsInSamples,ExpandedVCF-method}

\docType{methods}
\alias{variantsInSamples-methods}
\title{Identify variants observed in samples}
\description{
Identifies variants observed (uniquely) in at least one sample of a given
group.
}

\usage{
\S4method{variantsInSamples}{ExpandedVCF}(
    vcf, samples = 1:ncol(vcf), unique = FALSE)
}

\arguments{
    \item{vcf}{
\code{ExpandedVCF} object.

\code{metadata(vcf)[["TVTBparam"]]} must contain a
\code{\linkS4class{TVTBparam}} object.
}
    \item{samples}{
\code{integer}, \code{numeric} or \code{character} vector
indicating samples to consider in \code{VariantAnnotation::geno(vcf)}.
If not specified, all samples are considered.
}
    \item{unique}{
If \code{TRUE}, consider only variants \emph{unique} to the phenotype level
(\emph{i.e.} not seen in any other phenotype level).
}
}

\value{
An named \code{integer} vector of indices indicating the name and index of
variants that are (uniquely) observed in at least one non-reference genotype
in the given group of samples.
}

\section{Warning}{
A \code{warning} message is issued if genotypes are not fully defined in the
\code{TVTBparam}.
}

\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\linkS4class{VCF}}
and \code{\linkS4class{TVTBparam}}.
}

\examples{
# Example data ----

# VCF file
vcfFile <- system.file("extdata", "moderate.vcf", package = "TVTB")

# Phenotype file
phenoFile <- system.file("extdata", "moderate_pheno.txt", package = "TVTB")
phenotypes <- S4Vectors::DataFrame(
    read.table(file = phenoFile, header = TRUE, row.names = 1))

# TVTB parameters
tparam <- TVTBparam(Genotypes("0|0", c("0|1", "1|0"), "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(
    vcfFile, param = tparam, colData = phenotypes)
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)


# Example usage ----

variantsInSamples(
    vcf,
    which(SummarizedExperiment::colData(vcf)[,"super_pop"] == "EUR"))
}

\keyword{methods}

