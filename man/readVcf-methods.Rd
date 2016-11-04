\name{readVcf-methods}

\alias{readVcf,character,TVTBparam-method}
\alias{readVcf,TabixFile,TVTBparam-method}

\docType{methods}

\title{
Read VCF files
}

\description{
Read Variant Call Format (VCF) files, attaches the given \code{TVTBparam}
in the \code{metadata} slot of the resulting \code{VCF} object,
and attaches optional phenotype information in the \code{phenoData} slot.
}

\usage{
\S4method{readVcf}{character,TVTBparam}(
    file, genome, param, ..., colData = DataFrame(), autodetectGT = FALSE)
\S4method{readVcf}{TabixFile,TVTBparam}(
    file, genome, param, ..., colData = DataFrame(), autodetectGT = FALSE)
}

\arguments{
    \item{file, genome}{
See \code{\link{readVcf}}.
}
    \item{param}{
\code{TVTBparam} object that contains recurrent parameters.

The \code{vep} slot of \code{param} is checked for presence among the INFO
keys of the VCF file.
The \code{\linkS4class{TVTBparam}} object is coerced to
\code{\linkS4class{ScanVcfParam}} using the \code{ranges} slot only.
All \code{fixed}, \code{info}, and \code{geno} fields are imported
(see argument \code{colData} to declare samples to import).
}
    \item{...}{
Additional arguments, passed to methods.
}
    \item{colData}{
Phenotype information in a \code{\link{DataFrame}}.

If supplied, only samples identifiers present in \code{rownames(colData)}
are imported from the VCF file.
An error is thrown if any of the samples is absent from the VCF file.
}
    \item{autodetectGT}{
If \code{TRUE}, the method updates the genotypes definitions in the
\code{TVTBparam} object attached to the resulting \code{VCF} object after
guessing the codes that represent homozygote reference, heterozygote,
and homzoygote alternate genotypes.
}
}

\value{
\code{VCF} object.
See \code{?VCF} for complete details of the class structure.
}

\section{Warning}{
A \code{warning} message is issued if genotypes cannot be fully defined,
when \code{autodetectGT=TRUE}.
}

\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\link{readVcf,TabixFile,ScanVcfParam-method}},
and \code{\linkS4class{VCF}}.
}

\examples{
# Example data ----

# VCF file
vcfFile <- system.file("extdata", "moderate.vcf.gz", package = "TVTB")

# Phenotype file
phenoFile <- system.file("extdata", "moderate_pheno.txt", package = "TVTB")
phenotypes <- S4Vectors::DataFrame(read.table(phenoFile, TRUE, row.names = 1))

# TVTB parameters
tparam <- TVTBparam(Genotypes("0|0", c("0|1", "1|0"), "1|1"))


# Example usage ----

vcf <- readVcf(vcfFile, "b37", tparam, colData = phenotypes)
}
