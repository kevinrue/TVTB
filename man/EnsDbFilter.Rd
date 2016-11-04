\name{EnsDbFilter}
\alias{EnsDbFilter}

\title{
Dispatcher for construction of filter objects.
}
\description{
Create a filter object that can be used to retrieve specific elements from the
annotation database.
}
\usage{
EnsDbFilter(type, condition, value)
}

\arguments{
    \item{type}{
Type of filter. A prefix of the form "<type>Filter" on page
\code{\linkS4class{BasicFilter}}. See examples below for usage.
}
    \item{condition}{
The condition to be used in the comparison. See
\code{\linkS4class{BasicFilter}} for details.
}
    \item{value}{
The filter value. See \code{\linkS4class{BasicFilter}} for details.
}
}
\details{
Only \code{Genename} filter currently supported. Ideally, more filters and
combinations of filters will be supported in the future.
}
\value{
A \code{BasicFilter} that can be used to retrieve specific elements from the
annotation database.
}

\author{
Kevin Rue-Albrecht (adapted from Johannes Rainer)
}

\seealso{
\code{\linkS4class{BasicFilter}}
and \code{\link{genes,EnsDb-method}}
}

\examples{
    EnsDbFilter("Genename", "=", "SLC30A1")
}

