\name{tSVEParam-methods}

\alias{tSVEParam,list-method}
\alias{tSVEParam,missing-method}

\docType{methods}
\alias{tSVEParam-methods}
\title{Create tSVEParam objects}

\description{
Create a tSVEParam object
}

\usage{
\S4method{tSVEParam}{list}(
    genos, aaf, maf, bp)

\S4method{tSVEParam}{missing}(
    ref, het, alt, aaf, maf, bp)
}

\arguments{
    \item{genos}{
A (named) \code{list} of length three containing character vectors that declare homozygote reference, heterozygote, and homozygote alternate genotypes,
respectively.

If named, names set the individual INFO key suffix used to store the counts of
each set of genotypes. If unnamed, c("REF","HET","ALT") are used as default key suffixes.
}
    \item{ref}{
\code{character} vector of valid reference homozygote genotypes (\emph{e.g.}
    \code{c("0|0", "0/0")}).
}
    \item{het}{
\code{character} vector of valid heterozygote genotypes (\emph{e.g.}
    \code{c("0|1", "0/1")}).
}
    \item{alt}{
\code{character} vector of valid alternate homozygote genotypes (\emph{e.g.}
    \code{c("1|1", "1/1")}).
}
    \item{aaf}{
\code{character} vector of length 1 declaring the INFO key suffix used to
    store the alternate allele frequency (AAF).
}
    \item{maf}{
\code{character} vector of length 1 declaring the INFO key suffix used to
store the minor allele frequency (MAF).
}
    \item{bp}{
\code{BiocParallel} object.
}
}

\value{
A \code{tSVEParam} object that contains recurrent parameters.
}
\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\link{addFrequencies,ExpandedVCF,list,tSVEParam-method}},
\code{\link{tabulateCsqByPhenotype,ExpandedVCF,tSVEParam-method}},
and \code{\link{variantsInSamples,ExpandedVCF,numeric,tSVEParam-method}}.
variantsInSamples,ExpandedVCF,numeric,tSVEParam-method
}

\examples{
# Example usage ----
tparam <- tSVEParam(genos = list(
    REF = c("0|0"),
    HET = c("0|1","1|0"),
    ALT = c("1|1")))
}

\keyword{methods}

