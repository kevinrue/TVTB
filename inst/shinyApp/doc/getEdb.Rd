\name{getEdb}
\alias{getEdb}

\title{
Use EnsDb annotation packages
}
\description{
Based on the given \code{EnsDb} package name, load the library and return
the object.
}
\usage{
getEdb(x)
}

\arguments{
    \item{x}{
Character name of the \code{EnsDb} package.
}
}

\value{
\code{EnsDb} database object.
}

\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\linkS4class{EnsDb}}.
}
\examples{
if (require(EnsDb.Hsapiens.v75)){
    getEdb("EnsDb.Hsapiens.v75")
}
}


