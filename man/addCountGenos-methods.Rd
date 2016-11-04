\name{addCountGenos-methods}

\alias{addCountGenos}
\alias{addCountGenos,ExpandedVCF,missing-method}
\alias{addCountGenos,ExpandedVCF,numeric-method}
\alias{addCountGenos,ExpandedVCF,character-method}

\docType{methods}
\alias{addCountGenos-methods}
\title{Add count of genotypes to INFO field}
\description{
Adds the total occurences of a set of genotypes as an INFO
field for each variant.
All given genotypes are counted toward a single
total (\emph{e.g.} grand total of \code{c("0/0", "0|0")}),
while other genotypes are silently ignored.
}

\usage{
\S4method{addCountGenos}{ExpandedVCF,numeric}(
        vcf, genos, key, description, samples = 1:ncol(vcf), force = FALSE)

\S4method{addCountGenos}{ExpandedVCF,character}(
        vcf, genos, key, description, samples = 1:ncol(vcf), force = FALSE)
}

\arguments{
    \item{vcf}{
\code{ExpandedVCF} object.
}
    \item{genos}{
\code{character} vector of genotypes to count (toward a common unique total).
}
    \item{key}{
\code{character} name of the INFO field to create or update.
}
    \item{description}{
\code{character} description of the INFO field to create or update.
}
    \item{samples}{
\code{numeric} or \code{character} vector of indicating samples to consider
in \code{VariantAnnotation::geno(vcf)}.
}
    \item{force}{
If \code{TRUE}, the field header and data will be overwritten.

If \code{FALSE}, an error is thrown if the field already exists.
}
}

\value{
\code{ExpandedVCF} object including an additional \emph{INFO} field stating
    the count of genotypes.
}
\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\link{countGenos,ExpandedVCF,character-method}},
\code{\linkS4class{VCF}},
and \code{\link{geno,VCF-method}}
}

\examples{
# Example data ----

# VCF file
extdata <- file.path(system.file(package = "TVTB"), "extdata")
vcfFile <- file.path(extdata, "moderate.vcf")

# TVTB parameters
tparam <- TVTBparam(
    genos = list(
        REF = "0|0",
        HET = c("0|1", "1|0"),
        ALT = "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(file = vcfFile)
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)

# Example usage ----

## NOTE: use.names=FALSE is cleaner but not strictly necessary
addCountGenos(
    vcf = vcf, genos = unlist(het(tparam), use.names = FALSE),
    key = names(het(tparam)),
    description = "Number of heterozygous genotypes")
}

\keyword{methods}

