\name{VcfFilterList-class}
\docType{class}

% Class
\alias{class:VcfFilterList}
\alias{VcfFilterList-class}

% Constructors/Methods:
\alias{initialize,VcfFilterList-method}

\alias{VcfFilterList}
\alias{VcfFilterList-method}
\alias{VcfFilterList,ANY-method}

% Accessor methods:
\alias{filterRules}
\alias{filterRules,VcfFilterList-method}
\alias{filterRules<-}
\alias{filterRules<-,VcfFilterList,VcfFilterList-method}
\alias{filterRules<-,VcfFilterList,VcfFixedFilter-method}
\alias{filterRules<-,VcfFilterList,VcfInfoFilter-method}
\alias{filterRules<-,VcfFilterList,VcfVepFilter-method}
\alias{active}
\alias{active,VcfFilterList-method}
\alias{active<-}
\alias{active<-,VcfFilterList-method}


% Other methods / functions

\title{VcfFilterList class objects}

\description{
The \code{VcfFilterList} class stores a list of \code{VcfBasicFilter} objects.
}

\usage{
## Constructors

\S4method{VcfFilterList}{ANY}(..., active = rep(TRUE, length(list(...))))
}

\section{Accessors}{
In the following code snippets \code{x} is a \code{VcfFilterList} object.
\describe{
    \item{\code{filterRules(x)}, \code{filterRules(x) <- value}}{
    Gets or sets the name of the field to filter.
    }
    \item{\code{active(x)}, \code{active(x) <- value}}{
    Gets or sets the status of each filter.
    }
}
}

\arguments{
    \item{...}{
\code{VcfBasicFilter} objects.
}
    \item{active}{
A \code{logical} vector indicating active (\code{TRUE})
and inactive (\code{FALSE}) filters.
}
}

\value{
A \code{VcfFilterList} object.
}

\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\linkS4class{VcfBasicFilter}},
\code{\linkS4class{VCF}},
and \code{\link{parseCSQToGRanges,VCF-method}}.
}

\examples{
# Constructors ----

FF <- VcfFixedFilter(name = "FILTER", condition = "==", value = "PASS")
FF

IF <- VcfInfoFilter(name = "MAF", condition = "<=", value = 1e-6)
IF

VF <- VcfVepFilter(
    name = "IMPACT", condition = "\%in\%", value = c("MODERATE", "HIGH"))
VF

FL <- VcfFilterList(FF, IF, VF)
FL

# Accessors ----

## VcfBasicFilter objects stored in the VcfFilterList object
filterRules(FL)
active(FL)

# Setters ----
active(FL)[2] <- FALSE
FL
}
