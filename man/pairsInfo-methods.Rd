\name{pairsInfo-methods}

\alias{pairsInfo}
\alias{pairsInfo,VCF-method}

\docType{methods}
\alias{pairsInfo-methods}
\title{
Plot an INFO metric on a genomic axis.
}

\description{
Make a matrix of plots that display a metric calculated in levels of a given
phenotype,
and stored in columns of the \code{info} slot of a \code{VCF} object.
}

\usage{
\S4method{pairsInfo}{VCF}(vcf, metric, phenotype, ..., title = metric)
}

\arguments{
    \item{vcf}{
\code{VCF} object.
}
    \item{metric}{
Metric to plot on the Y axis.
All columns in the \code{info} slot of hte \code{vcf} object that match the
pattern \code{"phenotype_(.*)_metric"} are plotted in the \code{DataTrack}.
An error is thrown if no such column is found.
}
    \item{phenotype}{
Column in the \code{\link[Biobase]{phenoData}} slot of the \code{vcf} object.
Levels of this phenotype are plotted and contrasted in the \code{DataTrack}.
See argument \code{metric} for details.
}
    \item{...}{
Additional arguments, passed to the \code{\link[GGally]{ggpairs}} method.
}
    \item{title}{
Title for the graph, passed to the \code{ggpairs} method.
}
}

\value{
\code{gg} object returned by the \link[GGally]{ggpairs} method.
}

\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\link{ggpairs}},
\code{\link{addPhenoLevelFrequencies,ExpandedVCF-method}},
and \code{\linkS4class{VCF}}.
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
vcf <- addFrequencies(vcf, "super_pop")

# Example usage ----

pairsInfo(vcf, "MAF", "super_pop")
}
