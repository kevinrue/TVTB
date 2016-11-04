\name{preprocessVariants-methods}

\alias{preprocessVariants}

\alias{preprocessVariants,TabixFile,tSVEParam,DataFrame-method}
\alias{preprocessVariants,TabixFile,tSVEParam,data.frame-method}
\alias{preprocessVariants,TabixFile,tSVEParam,character-method}
\alias{preprocessVariants,TabixFile,tSVEParam,missing-method}

\alias{preprocessVariants,TabixFile,missing,DataFrame-method}
\alias{preprocessVariants,TabixFile,missing,data.frame-method}
\alias{preprocessVariants,TabixFile,missing,character-method}
\alias{preprocessVariants,TabixFile,missing,missing-method}

\alias{preprocessVariants,character,tSVEParam,DataFrame-method}
\alias{preprocessVariants,character,tSVEParam,data.frame-method}
\alias{preprocessVariants,character,tSVEParam,character-method}
\alias{preprocessVariants,character,tSVEParam,missing-method}

\alias{preprocessVariants,character,missing,DataFrame-method}
\alias{preprocessVariants,character,missing,data.frame-method}
\alias{preprocessVariants,character,missing,character-method}
\alias{preprocessVariants,character,missing,missing-method}

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
\S4method{preprocessVariants}{TabixFile,tSVEParam,DataFrame}(
    file, param, phenos, ...)

\S4method{preprocessVariants}{TabixFile,tSVEParam,data.frame}(
    file, param, phenos, ...)

\S4method{preprocessVariants}{TabixFile,tSVEParam,character}(
    file, param, phenos, ...)

\S4method{preprocessVariants}{TabixFile,tSVEParam,missing}(
    file, param, phenos = DataFrame(), ...)

\S4method{preprocessVariants}{TabixFile,missing,DataFrame}(
    file,
    param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
    phenos, ...)

\S4method{preprocessVariants}{TabixFile,missing,data.frame}(
    file,
    param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
    phenos, ...)

\S4method{preprocessVariants}{TabixFile,missing,character}(
    file,
    param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
    phenos, ...)

\S4method{preprocessVariants}{TabixFile,missing,missing}(
    file,
    param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
    phenos = DataFrame(), ...)

\S4method{preprocessVariants}{character,tSVEParam,DataFrame}(
    file, param, phenos, ...)

\S4method{preprocessVariants}{character,tSVEParam,data.frame}(
    file, param, phenos, ...)

\S4method{preprocessVariants}{character,tSVEParam,character}(
    file, param, phenos, ...)

\S4method{preprocessVariants}{character,tSVEParam,missing}(
    file, param, phenos = DataFrame(), ...)

\S4method{preprocessVariants}{character,missing,DataFrame}(
    file,
    param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
    phenos, ...)

\S4method{preprocessVariants}{character,missing,data.frame}(
    file,
    param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
    phenos, ...)

\S4method{preprocessVariants}{character,missing,data.frame}(
    file,
    param = tSVEParam(ref = "refNA", het = "hetNA", alt = "altNA"),
    phenos, ...)
}

\arguments{
    \item{file}{
File (\code{TabixFile} or \code{character}) from which to import variants.
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
vcfFile <- file.path(extdata, "moderate.vcf")
phenoFile <- file.path(extdata, "moderate_pheno.txt")

# Import phenotype data ----

phenotypes <- S4Vectors::DataFrame(read.table(
    file = phenoFile,
    header = TRUE,
    row.names = 1))

# Define tSVEParam

tparam <- tSVEParam(
    genos = list(
    REF = "0|0",
    HET = c("0|1", "1|0"),
    ALT = "1|1"))

# Example usage ----

preprocessVariants(
    file = vcfFile, param = tparam, phenos = phenotypes)
}

\keyword{methods}

