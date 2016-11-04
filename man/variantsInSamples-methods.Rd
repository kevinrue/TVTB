\name{variantsInSamples-methods}

\alias{variantsInSamples}
\alias{variantsInSamples,ExpandedVCF,numeric,TVTBparam-method}
\alias{variantsInSamples,ExpandedVCF,numeric,missing-method}
\alias{variantsInSamples,ExpandedVCF,character,TVTBparam-method}
\alias{variantsInSamples,ExpandedVCF,character,missing-method}

\docType{methods}
\alias{variantsInSamples-methods}
\title{Identify variants observed in samples}
\description{
Identifies variants observed (uniquely) in at least one sample of a given
group.
}

\usage{
\S4method{variantsInSamples}{ExpandedVCF,numeric,TVTBparam}(
    vcf, samples, param, ..., unique = FALSE)

\S4method{variantsInSamples}{ExpandedVCF,character,TVTBparam}(
    vcf, samples, param, ..., unique = FALSE)

\S4method{variantsInSamples}{ExpandedVCF,numeric,missing}(
    vcf, samples, alts, unique = FALSE)

\S4method{variantsInSamples}{ExpandedVCF,character,missing}(
    vcf, samples, alts, unique = FALSE)
}

\arguments{
    \item{vcf}{
\code{ExpandedVCF} object.
}
    \item{samples}{
\code{numeric} or \code{character} vector of indicating samples to consider
in \code{VariantAnnotation::geno(vcf)}.
}
    \item{param}{
\code{TVTBparam} object that contains recurrent parameters.
}
    \item{alts}{
\code{character} vector of alternate genotypes (both heterozygote and
homozygote).
}
    \item{...}{
Additional arguments used to override \code{TVTBparam} slots.
}
    \item{unique}{
If \code{TRUE}, consider only variants \emph{unique} to the phenotype level
(\emph{i.e.} not seen in any other phenotype level).
}
}

\value{
A numeric vector of indices indicating which variants are observed in at least
one non-reference genotype in the given group.
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
vcf <- VariantAnnotation::expand(vcf)

# Example usage ----

variantsInSamples(
    vcf = vcf,
    samples = which(SummarizedExperiment::colData(vcf)[,"super_pop"]=="EUR"),
    param = tparam)
}

\keyword{methods}

