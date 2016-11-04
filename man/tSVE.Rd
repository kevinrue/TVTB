\name{tSVE}
\alias{tSVE}

\title{
The Shiny Variant Explorer (tSVE) application
}
\description{
This function starts the interactive \code{tSVE} shiny web application that
    allows to interactively load and visualise genetic variants and their
    consequences predicted by the Ensembl VEP.
}
\usage{
    tSVE(...)
}

\arguments{
    \item{\dots}{
Additional arguments passed to the \code{\link[shiny]{runApp}} function
from the \code{shiny} package.
}
}

\references{
Interface to \code{EnsDb} adapted from \code{ensembldb}.
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
