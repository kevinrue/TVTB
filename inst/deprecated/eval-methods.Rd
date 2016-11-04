\name{eval-methods}

\alias{eval}

\alias{eval,VcfFixedRules,ExpandedVCF-method}
\alias{eval,VcfInfoRules,ExpandedVCF-method}
\alias{eval,VcfVepRules,ExpandedVCF-method}
\alias{eval,VcfFilterRules,ExpandedVCF-method}

\docType{methods}
\alias{eval-methods}
\title{Count occurences of genotypes}
\description{
\code{VcfFixedRules} and \code{VcfFixedRules} are collections of filter rules
applicable to their corresponding \code{slot} in an \code{ExpandedVCF} object.
\code{VcfFixedRules} filter rules are applicable to the Ensembl VEP predictions
stored in a given INFO field of an \code{ExpandedVCF} object.
\code{VcfFilterRules} are collection of the above filter rules, each marked
to their respective target slot/object.
}

\usage{
\S4method{eval}{VcfFixedRules,ExpandedVCF}(expr, envir)

\S4method{eval}{VcfInfoRules,ExpandedVCF}(expr, envir)

\S4method{eval}{VcfVepRules,ExpandedVCF}(expr, envir)

\S4method{eval}{VcfFilterRules,ExpandedVCF}(expr, envir)
}

\arguments{
    \item{expr}{
\code{VcfFixedRules}, \code{VcfFixedRules}, \code{VcfFixedRules},
or\code{VcfFilterRules} object.
}
    \item{envir}{
\code{ExpandedVCF} object.
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
# Constructors ----

fixedR <- VcfFixedRules(exprs = list(
    pass = expression(FILTER == "PASS"),
    qual = expression(QUAL > 20)
    ))
fixedR

infoR <- VcfInfoRules(exprs = list(
    common = expression(MAF > 0.1), # minor allele frequency
    alt = expression(ALT > 0) # count of alternative homozygotes
    ))
infoR

vepR <- VcfVepRules(exprs = list(
    missense = expression(Consequence \%in\% c("missense_variant")),
    CADD = expression(CADD_PHRED > 15)
    ))
vepR

vcfRules <- VcfFilterRules(fixedR, infoR, vepR)
vcfRules

# Example data ----

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
vcf <- VariantAnnotation::expand(vcf)
vcf <- addOverallFrequencies(vcf = vcf, param = tparam)

# Applying filters to ExpandedVCF objects ----

## Evaluate filters
S4Vectors::eval(expr = fixedR, envir = vcf)
S4Vectors::eval(expr = infoR, envir = vcf)
S4Vectors::eval(expr = vepR, envir = vcf)
S4Vectors::eval(expr = vcfRules, envir = vcf)

## Evaluate filters separately
as.data.frame(S4Vectors::evalSeparately(expr = vepR, envir = vcf))
as.data.frame(S4Vectors::evalSeparately(expr = vcfRules, envir = vcf))

## Subset ExpandedVCF by filters
S4Vectors::subsetByFilter(x = vcf, filter = vepR)
S4Vectors::subsetByFilter(x = vcf, filter = vcfRules)
}

\keyword{methods}

