\name{tSVE}
\alias{tSVE}

\title{
The Shiny Variant Explorer (tSVE) web-application
}
\description{
\emph{Currently unsupported --- Package undergoing major updates.}

This function starts the interactive \code{tSVE} shiny web-application that
allows to interactively load and visualise genetic variants and their
Ensembl Variant Effect Predictor (VEP) predictions using the package methods.
All arguments after the \dots set default values for the application
(\emph{e.g.} widgets).
}
\usage{
tSVE(
    ...,
    refGT = "0|0",
    hetGT = c("0|1", "1|2", "0|2", "1|0", "2|1", "2|0"),
    altGT = c("1|1", "2|2"),
    vepKey = "CSQ",
    refSuffix = "REF", hetSuffix = "HET", altSuffix = "ALT",
    aafSuffix = "AAF", mafSuffix = "MAF",
    genoHeatmap.height = "500px",
    options.width = 120,
    autodetectGTimport = FALSE
)
}

\arguments{
    \item{\dots}{
Additional arguments passed to the \code{\link[shiny]{runApp}} function
from the \code{shiny} package.
}
    \item{refGT}{
Default homozygote reference genotypes.
}
    \item{hetGT}{
Default heterozygote genotypes.
}
    \item{altGT}{
Default homozygote alternate genotypes.
}
    \item{vepKey}{
Default INFO key for the VEP prediction field.
}
    \item{refSuffix}{
Default INFO key suffix used to store the data for homozygote reference
genotypes.
}
    \item{hetSuffix}{
Default INFO key suffix used to store the data for heterozygote genotypes.
}
    \item{altSuffix}{
Default INFO key suffix used to store the data for homozygote alternate
genotypes.
}
    \item{aafSuffix}{
Default INFO key suffix used to store the data for alternate allele frequency.
}
    \item{mafSuffix}{
Default INFO key suffix used to store the data for minor allele frequency.
}
    \item{genoHeatmap.height}{
Default height (in pixels) of the heatmap that represents the genotype of
each variant in each sample.
}
    \item{options.width}{
Sets \code{options("width")}.
}
    \item{autodetectGTimport}{
Default checkbox value. If \code{FALSE}, genotypes (ref, het, alt) are
taken \emph{as is} from the \emph{Advanced settings} panel. If \code{TRUE},
genotypes selected in the \emph{Advanced settings} panel are updated using
the \code{autodetectGenotypes} method, immediately after variants are imported.
}

}

\references{
Interface to \code{EnsDb} adapted from the \code{ensembldb} package.
}
\value{
Not applicable (yet).
}
\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\link{runEnsDbApp}}.
}
\examples{
if (interactive()){
    runEnsDbApp()
}
}
