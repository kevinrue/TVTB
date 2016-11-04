\name{preprocessVariants}

\title{
Import and pre-process variants from VCF file(s).
}
\description{
\code{mafFilter()} returns a \code{ExpandedVCF} object retaining variants
    within a range of allele frequencies.
}

\usage{
mafFilter(vcf, MAF)
}

\arguments{
    \item{vcf}{
\code{ExpandedVCF} object.
}
    \item{MAF}{
Numeric vector of length 2 specyfing minimum and maximum minor allele frequency
    to filter imported variants, respectively.
    \emph{Values must be in range [0,0.5]}.
}
}

\value{
\describe{
\item{mafFilter}{\code{ExpandedVCF} object including variants within the
    given range.}
}
}

\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\link{addFrequencies}},
\code{\link{countGenos}},
\code{\link{chr2file}},
\code{\link{readVcf}},
\code{\linkS4class{VCF}},
\code{\link{getEdb}},
\code{\linkS4class{TabixFile}},
\code{\linkS4class{GRanges}},
\code{\linkS4class{Seqinfo}},
and \code{\linkS4class{BiocParallelParam}}.
}

\examples{
# Example data ----

extdata <- file.path(system.file(package = "tSVE"), "extdata")
phenoFile <- file.path(extdata, "integrated_samples.txt")
regions <- GenomicRanges::GRanges(
    seqnames = "15",
    ranges = IRanges::IRanges(start = 48420E3, end = 48421E3))

vcfFile <- file.path(extdata, "chr15.phase3_integrated.vcf.gz")

# Import phenotype data ----

phenotypes <- S4Vectors::DataFrame(read.table(
    file = phenoFile,
    header = TRUE,
    row.names = 1))

# Reference Tabix file ----

tabixVcf <- Rsamtools::TabixFile(vcfFile)

# Example usage ----

vcf <- preprocessVariants(
    file = tabixVcf, regions = regions,
    ref = c("0|0"), het = c("0|1", "1|0"), alt = c("1|1"),
    phenos = phenotypes)
}

