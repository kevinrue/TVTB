\name{countGenos-methods}

\alias{countGenos}

\alias{countGenos,ExpandedVCF,character-method}
\alias{countGenos,matrix,character-method}

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
\S4method{countGenos}{ExpandedVCF,character}(
    x, genos, pheno = NULL, level = NULL)

\S4method{countGenos}{matrix,character}(
    x, genos)
}

\arguments{
    \item{x}{
\code{ExpandedVCF} or \code{matrix} of genotypes
(\emph{e.g.} \code{geno(vcf)[["GT"]]}).
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

## NOTE: use.names=FALSE is cleaner but not strictly necessary
countGenos(
    x = vcf,
    genos = unlist(het(tparam), use.names = FALSE),
    pheno = "super_pop",
    level = "AFR")
}

\keyword{methods}

