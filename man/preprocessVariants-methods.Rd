\name{preprocessVariants-methods}

\alias{preprocessVariants}

\alias{preprocessVariants,TabixFile,GRanges,tSVEParam,DataFrame-method}
\alias{preprocessVariants,TabixFile,GRanges,tSVEParam,data.frame-method}
\alias{preprocessVariants,TabixFile,GRanges,tSVEParam,character-method}
\alias{preprocessVariants,TabixFile,GRanges,tSVEParam,missing-method}

\alias{preprocessVariants,TabixFile,GRanges,missing,DataFrame-method}
\alias{preprocessVariants,TabixFile,GRanges,missing,data.frame-method}
\alias{preprocessVariants,TabixFile,GRanges,missing,character-method}
\alias{preprocessVariants,TabixFile,GRanges,missing,missing-method}

\alias{preprocessVariants,character,GRanges,tSVEParam,DataFrame-method}
\alias{preprocessVariants,character,GRanges,tSVEParam,data.frame-method}
\alias{preprocessVariants,character,GRanges,tSVEParam,character-method}
\alias{preprocessVariants,character,GRanges,tSVEParam,missing-method}

\alias{preprocessVariants,character,GRanges,missing,DataFrame-method}
\alias{preprocessVariants,character,GRanges,missing,data.frame-method}
\alias{preprocessVariants,character,GRanges,missing,character-method}
\alias{preprocessVariants,character,GRanges,missing,missing-method}

\docType{methods}
\alias{preprocessVariants-methods}
\title{Pre-process genetic variants from VCF file}
\description{
Imports variants from a VCF file
within given genomic regions, along with their Ensembl VEP predictions,
to initialise a minimal \code{ExpandedVCF} object required for subsequent
analyses.
}

\usage{
\S4method{preprocessVariants}{TabixFile,GRanges,tSVEParam,DataFrame}(
    file, regions, param, phenos, ...)

\S4method{preprocessVariants}{TabixFile,GRanges,tSVEParam,data.frame}(
    file, regions, param, phenos, ...)

\S4method{preprocessVariants}{TabixFile,GRanges,tSVEParam,character}(
    file, regions, param, phenos, ...)

\S4method{preprocessVariants}{TabixFile,GRanges,tSVEParam,missing}(
    file, regions, param, phenos = DataFrame(), ...)

\S4method{preprocessVariants}{TabixFile,GRanges,missing,DataFrame}(
    file, regions,
    param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
    phenos, ...)

\S4method{preprocessVariants}{TabixFile,GRanges,missing,data.frame}(
    file, regions,
    param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
    phenos, ...)

\S4method{preprocessVariants}{TabixFile,GRanges,missing,character}(
    file, regions,
    param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
    phenos, ...)

\S4method{preprocessVariants}{TabixFile,GRanges,missing,missing}(
    file, regions,
    param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
    phenos = DataFrame(), ...)

\S4method{preprocessVariants}{character,GRanges,tSVEParam,DataFrame}(
    file, regions, param, phenos, ...)

\S4method{preprocessVariants}{character,GRanges,tSVEParam,data.frame}(
    file, regions, param, phenos, ...)

\S4method{preprocessVariants}{character,GRanges,tSVEParam,character}(
    file, regions, param, phenos, ...)

\S4method{preprocessVariants}{character,GRanges,tSVEParam,missing}(
    file, regions, param, phenos = DataFrame(), ...)

\S4method{preprocessVariants}{character,GRanges,missing,DataFrame}(
    file, regions,
    param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
    phenos, ...)

\S4method{preprocessVariants}{character,GRanges,missing,data.frame}(
    file, regions,
    param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
    phenos, ...)

\S4method{preprocessVariants}{character,GRanges,missing,data.frame}(
    file, regions,
    param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
    phenos, ...)

\S4method{preprocessVariants}{character,GRanges,missing,missing}(
    file, regions,
    param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
    phenos = DataFrame(), ...)
}

\arguments{
    \item{file}{
File (\code{TabixFile} or \code{character}) from which to import variants.
}
    \item{regions}{
\code{GRanges} object specifying genomic regions from which to import variants.
}
    \item{param}{
\code{tSVEParam} object that contains recurrent parameters.
}
    \item{phenos}{
\code{DataFrame}, \code{data.frame}, or file of phenotypes.

In all cases, \code{rownames} must provide sample names,
and \code{colnames} must provide phenotype names.
}
    \item{...}{
Additional arguments passed to \code{tSVEParam} constructor and other methods
(\emph{e.g.} argument \code{genome} for method \code{readVcf}).
}
}

\value{
\code{ExpandedVCF} object pre-processed as follows:
\enumerate{
    \item{imported by \code{\link{readVcf}},}
    \item{expanded to bi-allelic records,}
    \item{combined with optional \code{phenos} information;
    \emph{Highly recommended for subsequent analyses}}
}
}

\details{
If \code{vcfTabix} contains multi-allelic records, \code{row.names} of the
\code{ExpandedVCF} object are disambiguated by appending the ALT allele to
the original \code{row.names}.
}

\note{
Currently, \code{readVCF} is used to import only:
\describe{
    \item{fixed}{ALT field}
    \item{info}{\code{vep} field}
    \item{geno}{GT field}
}
}

\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\linkS4class{TabixFile}},
\code{\linkS4class{GRanges}},
\code{\link{readVcf}},
\code{\linkS4class{VCF}},
\code{\link{fixed,VCF-method}},
\code{\link{info,VCF-method}},
\code{\link{geno,VCF-method}}
}

\examples{
library(VariantAnnotation)

# Example data ----

extdata <- file.path(system.file(package = "tSVE"), "extdata")
vcfFile <- file.path(extdata, "chr15.phase3_integrated.vcf.gz")
phenoFile <- file.path(extdata, "integrated_samples.txt")
regions <- GenomicRanges::GRanges(
    seqnames = "15",
    ranges = IRanges::IRanges(start = 48420E3, end = 48421E3))

# Import phenotype data ----

phenotypes <- S4Vectors::DataFrame(read.table(
    file = phenoFile,
    header = TRUE,
    row.names = 1))

# Define VCF and index files ----

tf <- Rsamtools::TabixFile(vcfFile)

# Define tSVEParam

tParam <- tSVEParam(
    ref = "0|0", het = c("0|1", "1|0"), alt = "1|1",
    vep = "CSQ")

# Example usage ----

preprocessVariants(
    file = tf, regions = regions, param = tParam,
    phenos = phenotypes)
}

\keyword{methods}

