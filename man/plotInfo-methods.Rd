\name{plotInfo-methods}

\alias{plotInfo}
\alias{plotInfo,VCF-method}

\docType{methods}
\alias{plotInfo-methods}
\title{
Plot an INFO metric on a genomic axis.
}

\description{
Plot, on a genomic axis, a metric calculated in levels of a given phenotype,
and stored in columns of the \code{info} slot of a \code{VCF} object.
}

\usage{
\S4method{plotInfo}{VCF}(
    vcf, metric, range, annotation, phenotype, type = c("p", "heatmap"))
}

\arguments{
    \item{vcf}{
\code{ExpandedVCF} object.
}
    \item{metric}{
Metric to plot on the Y axis.
All columns in the \code{info} slot of hte \code{vcf} object that match the
pattern \code{"phenotype_(.*)_metric"} are plotted in the \code{DataTrack}.
An error is thrown if no such column is found.
}
    \item{range}{
A \code{GRanges} of length one that defines the region to visualise. All
variants in the \code{vcf} object overlapping this region are plotted.
}
    \item{annotation}{
An \code{\linkS4class{EnsDb}} annotation package from which to fetch gene
annotations.
\emph{\code{\link[GenomicFeatures]{TxDb}} packages may be supported in the
future}.
}
    \item{phenotype}{
Column in the \code{\link[Biobase]{phenoData}} slot of the \code{vcf} object.
Levels of this phenotype are plotted and contrasted in the \code{DataTrack}.
See argument \code{metric} for details.
}
    \item{type}{
Plotting type(s), as listed in \code{\linkS4class{DataTrack}}.
}
}

\value{
\code{list} returned by the \code{plotTracks} method.
}

\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\link{plotTracks}},
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

if (requireNamespace("EnsDb.Hsapiens.v75")){
    plotInfo(
        vcf, "MAF",
        range(GenomicRanges::granges(vcf)),
        EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
        "super_pop"
    )
}


}
