\name{countGenos-methods}

\alias{countGenos}

\alias{countGenos,ExpandedVCF-method}

\docType{methods}
\alias{countGenos-methods}
\title{Count occurences of genotypes}
\description{
Counts the total occurences of a set of genotypes by row
in a matrix of genotype. All given genotypes are counted toward a single
total (\emph{e.g.} grand total of \code{c("0/0", "0|0")}),
while other genotypes are silently ignored.
}

\usage{
\S4method{countGenos}{ExpandedVCF}(
    x, genos, pheno = NULL, level = NULL)
}

\arguments{
    \item{x}{
\code{ExpandedVCF} object.
}
    \item{genos}{
\code{character} vector of genotypes to count (toward a common unique total).
}
    \item{pheno}{
If \code{x} is an \code{ExpandedVCF} object, phenotype in
\code{colnames(colData(x))}.
}
    \item{level}{
If \code{x} is an \code{ExpandedVCF} object, phenotype level in
\code{colData(x)[,pheno]}.
}
}

\value{
An \code{integer} vector representing the aggregated count of the given
genotypes in each row.
}
\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\linkS4class{VCF}}
}

\examples{
# Example data ----

# VCF file
vcfFile <- system.file("extdata", "moderate.vcf", package = "TVTB")

# Phenotype file
phenoFile <- system.file("extdata", "moderate_pheno.txt", package = "TVTB")
phenotypes <- S4Vectors::DataFrame(read.table(phenoFile, TRUE, row.names = 1))

# TVTB parameters
tparam <- TVTBparam(Genotypes("0|0", c("0|1", "1|0"), "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(
    vcfFile, param = tparam, colData = phenotypes)
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)


# Example usage ----

vcf <- countGenos(vcf, het(tparam), "super_pop", "AFR")
}

\keyword{methods}

