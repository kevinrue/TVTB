\name{subsetVcf-methods}

\alias{subsetVcf}
\alias{subsetVcf,ExpandedVCF,VcfFixedFilter,missing-method}
\alias{subsetVcf,ExpandedVCF,VcfInfoFilter,missing-method}
\alias{subsetVcf,ExpandedVCF,VcfVepFilter,TVTBparam-method}
\alias{subsetVcf,ExpandedVCF,VcfVepFilter,missing-method}
\alias{subsetVcf,ExpandedVCF,VcfFilterList,TVTBparam-method}
\alias{subsetVcf,ExpandedVCF,VcfFilterList,missing-method}

\docType{methods}
\alias{subsetVcf-methods}
\title{
Subset VCF objects
}

\description{
Subset VCF objects using filters based on any of:
\itemize{
\item{\code{fixed} fields}
\item{\code{info} fields}
\item{Ensembl VEP fields \emph{(stored in a given INFO field)}}
}


}

\usage{
\S4method{subsetVcf}{ExpandedVCF,VcfFixedFilter,missing}(
    x, filter)

\S4method{subsetVcf}{ExpandedVCF,VcfInfoFilter,missing}(
    x, filter)

\S4method{subsetVcf}{ExpandedVCF,VcfVepFilter,TVTBparam}(
    x, filter, param, ..., vep = FALSE)

\S4method{subsetVcf}{ExpandedVCF,VcfVepFilter,missing}(
    x, filter,
    param = TVTBparam(ref = "refNA", het = "hetNA", alt = "altNA"),
    ..., vep = FALSE)
}

\arguments{
    \item{x}{
\code{ExpandedVCF} object.
}
    \item{filter}{
\code{VcfBasicFilter} object.
}
    \item{param}{
\code{TVTBparam} object that contains recurrent parameters.
In this case, the INFO key that contains Ensembl VEP predictions.
}
    \item{...}{
Additional arguments used to override \code{TVTBparam} slots.
}
    \item{vep}{
\strong{Not implemented yet}
If \code{TRUE}, Ensembl VEP predictions are also subsetted to those
passing the filters.

If \code{FALSE}, all Ensembl VEP predictions of the variants passing the
filters are kept.
}
}

\value{
\code{ExpandedVCF} object subsetted to variants passing the filters.
}

\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\linkS4class{VcfBasicFilter}}
and \code{\linkS4class{VCF}}.
}

\examples{
# VCF file
extdata <- file.path(system.file(package = "TVTB"), "extdata")
vcfFile <- file.path(extdata, "moderate.vcf")

# Parameters
tparam <- TVTBparam(
    genos = list(
        REF = "0|0",
        HET = c("0|1", "1|0"),
        ALT = "1|1"))

vcf <- VariantAnnotation::readVcf(file = vcfFile)
vcf <- VariantAnnotation::expand(x = vcf)
vcf <- addOverallFrequencies(vcf = vcf, param = tparam)

# Example usage ----

fixedFilter <- VcfFixedFilter(
    name = "FILTER", condition = "==", value = "PASS")
subsetVcf(x = vcf, filter = fixedFilter)

infoFilter <- VcfInfoFilter(name = "MAF", condition = ">", value = 0.1)
subsetVcf(x = vcf, filter = infoFilter)

vepFilter <- VcfVepFilter(
    name = "Consequence", condition = "\%in\%", value = c("missense_variant"))
subsetVcf(x = vcf, filter = vepFilter, param = tparam, vep = FALSE)

listFilter <- VcfFilterList(
    exprs = list(fixedFilter, infoFilter, vepFilter))
subsetVcf(x = vcf, filter = listFilter, param = tparam, vep = FALSE)
}
