\name{altFreq}
\alias{altFreq}

\title{
Calculate alternate allele frequency (AAF)
}
\description{
Calculate alternate allele frequency from matrix of genotypes and vectors of
homozygote and heterozygote genotypes.
}
\usage{
altFreq(GT, ref, het, alt)
}

\arguments{
    \item{GT}{
Matrix of genotypes (\emph{e.g.} \code{geno(vcf)[["GT"]]}).
}
    \item{ref}{
Vector of reference homozygote genotypes (\emph{e.g.} \code{c("0/0")}).
}
    \item{het}{
Vector of heterozygote genotypes (\emph{e.g.} \code{"0/1"})
}
    \item{alt}{
Vector of alternate homozygote genotypes (\emph{e.g.} \code{c("1/1")})
}
}

\value{
Numeric vector of alternate allele frequencies.
}

\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\linkS4class{VCF}}.
}
\examples{
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

svp <- VariantAnnotation::ScanVcfParam(
    geno = "GT",
    info = "CSQ",
    samples = rownames(phenotypes),
    which = bedRegions)
tf <- Rsamtools::TabixFile(vcfFile)
vcf <- VariantAnnotation::readVcf(file = tf, param = svp)

# Example usage ----

altFreq(
    GT = VariantAnnotation::geno(vcf)[["GT"]],
    ref = c("0|0"),
    het = c("0|1"),
    alt = c("1|1"))

}

