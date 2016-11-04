\name{autodetectGenotypes-methods}

\alias{autodetectGenotypes}
\alias{autodetectGenotypes,VCF-method}

\docType{methods}
\alias{autodetectGenotypes-methods}
\title{Define genotypes in the TVTBparam metadata slot}
\description{
This method attempts to auto-detect genotypes (\emph{i.e.}
homozygote reference, heterozygote, and homozygote alternate)
in a \code{VCF} object, and sets or creates a \code{TVTBparam}
object accordingly, in the \code{metadata} slot.
}

\usage{
\S4method{autodetectGenotypes}{VCF}(vcf)
}

\arguments{
    \item{vcf}{
\code{VCF} object.
}
}

\value{
\code{VCF} object including a new or updated \code{TVTBparam} object
in \code{metadata(vcf)[["TVTBparam"]]} .
}

\section{Warning}{
A \code{warning} message is issued if genotypes cannot be fully defined.
}

\author{
Kevin Rue-Albrecht
}

\examples{
# Example data ----

# VCF file
vcfFile <- system.file("extdata", "moderate.vcf", package = "TVTB")

# TVTB parameters
tparam <- TVTBparam()

# Pre-process variants
vcf <- VariantAnnotation::readVcf(vcfFile, param = tparam) # warning expected
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)


# Example usage ----

vcf <- autodetectGenotypes(vcf)
}

\keyword{methods}

