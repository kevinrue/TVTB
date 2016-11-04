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
# Separate multi-allelic records into bi-allelic records
eVcf <- expand(x = vcf, row.names = TRUE)
# Disambiguate row.names from multi-allelic records
rownames(eVcf) <- paste(rownames(eVcf), mcols(eVcf)[,"ALT"], sep = "_")

# Example usage ----
addCountGenos(
    vcf = eVcf, genos = c("0|1", "1|0"),
    key = "NHET",
    description = "Number of heterozygous genotypes")
}

\keyword{methods}

