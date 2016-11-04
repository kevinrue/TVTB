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
library(VariantAnnotation)

# Example data ----

extdata <- file.path(system.file(package = "tSVE"), "extdata")
vcfFile <- file.path(extdata, "chr15.phase3_integrated.vcf.gz")
phenoFile <- file.path(extdata, "integrated_samples.txt")
bedRegions <- GenomicRanges::GRanges(
    seqnames = "15",
    ranges = IRanges::IRanges(start = 48420E3, end = 48421E3))

# Import phenotype data ----

phenotypes <- S4Vectors::DataFrame(read.table(
    file = phenoFile,
    header = TRUE,
    row.names = 1))

# Import variants ----

svp <- ScanVcfParam(
    geno = "GT",
    info = "CSQ",
    samples = rownames(phenotypes),
    which = bedRegions)
tf <- Rsamtools::TabixFile(vcfFile)
vcf <- readVcf(file = tf, param = svp)
colData(vcf) <- phenotypes

# Separate multi-allelic records into bi-allelic records
eVcf <- expand(x = vcf, row.names = TRUE)
# Disambiguate row.names from multi-allelic records
rownames(eVcf) <- paste(rownames(eVcf), mcols(eVcf)[,"ALT"], sep = "_")

# Example usage ----

countGenos(x = eVcf, genos = c("0|1"), pheno = "super_pop", level = "AFR")
}

\keyword{methods}

