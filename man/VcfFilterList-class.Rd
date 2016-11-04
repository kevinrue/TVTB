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

% Other methods / functions

\title{VcfFilterList class objects}

\description{
The \code{VcfFilterList} class stores a list of \code{VcfBasicFilter} objects.
}

\usage{
## Constructors

\S4method{VcfFilterList}{ANY}(...)
}

\section{Accessors}{
In the following code snippets \code{x} is a \code{VcfFilterList} object.
\describe{
    \item{\code{filterRules(x)}, \code{filterRules(x) <- value}}{
    Gets or sets the name of the field to filter.
    }
}
}

\arguments{
    \item{...}{
\code{VcfBasicFilter} objects.
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

vff <- VcfFixedFilter(name = "FILTER", condition = "==", value = "PASS")
vff

vif <- VcfInfoFilter(name = "MAF", condition = "<=", value = 1e-6)
vif

vvf <- VcfVepFilter(
    name = "IMPACT", condition = "\%in\%", value = c("MODERATE", "HIGH"))
vvf

vlist <- VcfFilterList(vff, vif, vvf)
vlist

# Accessors ----

## VcfBasicFilter objects stored in the VcfFilterList object
filterRules(vlist)
}
