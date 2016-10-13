\name{addCountGenos-methods}

\alias{addCountGenos}
\alias{addCountGenos,ExpandedVCF-method}

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
\S4method{addCountGenos}{ExpandedVCF}(
        vcf, genos, key, description,
        samples = 1:ncol(vcf), force = FALSE)
}

\arguments{
    \item{vcf}{
\code{ExpandedVCF} object.
}
    \item{genos}{
\code{character} vector of genotypes to count (toward a common unique total).
}
    \item{key}{
Name of the INFO field to create or update
(\code{character} vector of length 1). See \emph{Details} below.
}
    \item{description}{
\code{character} description of the INFO field to create or overwrite
(\code{character} vector of length 1).
}
    \item{samples}{
\code{integer}, \code{numeric} or \code{character} vector
indicating samples to consider in \code{VariantAnnotation::geno(vcf)}.
If not specified, all samples are considered.
}
    \item{force}{
If \code{TRUE}, the field header and data will be overwritten if present;
If \code{FALSE}, an error is thrown if the field already exists.
}
}

\value{
\code{ExpandedVCF} object including an additional INFO field stating
    the count of genotypes.
}
\author{
Kevin Rue-Albrecht
}

\details{
In all cases, the new INFO field is inserted after the last existing field.
In other words, overwriting an existing INFO field is achieved by dropping
it from the data and header of the \code{info} slot,
and subsequently inserting the new data after the last remaining INFO field.
}

\seealso{
\code{\link{countGenos,ExpandedVCF-method}}
and \code{\link{geno,VCF-method}}
}

\examples{
# Example data ----

# VCF file
vcfFile <- system.file("extdata", "moderate.vcf", package = "TVTB")

# TVTB parameters
tparam <- TVTBparam(Genotypes(ref = "0|0", het = c("0|1", "1|0"), alt = "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(vcfFile, param = tparam)
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)


# Example usage ----

vcf <- addCountGenos(
    vcf, het(tparam),
    suffix(tparam)["het"],
    "Number of heterozygous genotypes")
}

\keyword{methods}

