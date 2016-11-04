\name{dropInfo-methods}

\alias{dropInfo}

\alias{dropInfo,VCF-method}

\docType{methods}
\alias{dropInfo-methods}

\title{Remove INFO keys from VCF objects}

\description{
Given a \code{character} vector of INFO keys, removes either the associated
header, data, or both from a \code{VCF} object.
If no INFO key is given (the default), all INFO keys are checked and removed
from the given slot if they do not have a matching entry in the other slot.
}

\usage{
\S4method{dropInfo}{VCF}(
    vcf, key = NULL, slot = "both")
}

\arguments{
    \item{vcf}{
\code{VCF} object.
}
    \item{key}{
\code{character} vector of INFO keys to remove.

If \code{NULL} (the default), all keys are checked,
and removed from the given \code{slot} if they do not have a matching entry
in the other slot.
}
    \item{slot}{
Should the INFO keys be removed from the "header", the "data",
or \code{"both"} (the default)?
}
}

\value{
An \code{integer} vector representing the aggregated count of the given
genotypes in each row.
}
\author{
Kevin Rue-Albrecht
}

\note{
In the future, \code{x} should also support genotype quality (GQ) to consider
only genotypes above a given quality cut-off.
}

\seealso{
\code{\linkS4class{VCF}}
}

\examples{
# Example data ----

# VCF file
vcfFile <- system.file("extdata", "moderate.vcf", package = "TVTB")

# TVTB parameters
tparam <- TVTBparam(Genotypes("0|0", c("0|1", "1|0"), "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(
    vcfFile, param = tparam)
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)


# Example usage ----

dropInfo(vcf)
dropInfo(vcf, "CSQ")
}

\keyword{methods}

