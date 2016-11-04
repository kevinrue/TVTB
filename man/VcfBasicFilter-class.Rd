\name{VcfBasicFilter-class}
\docType{class}

% Class
\alias{class:VcfBasicFilter}
\alias{VcfBasicFilter-class}
\alias{class:VcfFixedFilter}
\alias{VcfFixedFilter-class}
\alias{class:VcfInfoFilter}
\alias{VcfInfoFilter-class}
\alias{class:VcfVepFilter}
\alias{VcfVepFilter-class}

% Constructors/Methods:
\alias{initialize,VcfFixedFilter-method}
\alias{initialize,VcfInfoFilter-method}
\alias{initialize,VcfVepFilter-method}

\alias{VcfFixedFilter}
\alias{VcfFixedFilter-methods}
\alias{VcfFixedFilter,character,character-method}
\alias{VcfFixedFilter,character,character,ANY-method}

\alias{VcfInfoFilter}
\alias{VcfInfoFilter-methods}
\alias{VcfInfoFilter,character,character-method}
\alias{VcfInfoFilter,character,character,ANY-method}

\alias{VcfVepFilter}
\alias{VcfVepFilter-methods}
\alias{VcfVepFilter,character,character-method}
\alias{VcfVepFilter,character,character,ANY-method}

% Accessor methods:
\alias{name}
\alias{name<-}
\alias{condition}
\alias{condition<-}
\alias{value}
\alias{value<-}
\alias{filterType}

\alias{name,VcfFixedFilter-method}
\alias{name<-,VcfFixedFilter,character-method}
\alias{condition,VcfFixedFilter-method}
\alias{condition<-,VcfFixedFilter,character-method}
\alias{value,VcfFixedFilter-method}
\alias{value<-,VcfFixedFilter-method}
\alias{filterType,VcfFixedFilter-method}

\alias{name,VcfInfoFilter-method}
\alias{name<-,VcfInfoFilter,character-method}
\alias{condition,VcfInfoFilter-method}
\alias{condition<-,VcfInfoFilter,character-method}
\alias{value,VcfInfoFilter-method}
\alias{value<-,VcfInfoFilter-method}
\alias{filterType,VcfInfoFilter-method}

\alias{name,VcfVepFilter-method}
\alias{name<-,VcfVepFilter,character-method}
\alias{condition,VcfVepFilter-method}
\alias{condition<-,VcfVepFilter,character-method}
\alias{value,VcfVepFilter-method}
\alias{value<-,VcfVepFilter-method}
\alias{filterType,VcfVepFilter-method}

% Other methods / functions
\alias{as.character}
\alias{as.character,VcfFixedFilter-method}
\alias{as.character,VcfInfoFilter-method}
\alias{as.character,VcfVepFilter-method}

\title{VcfBasicFilter class objects}

\description{
The \code{VcfBasicFilter} class stores the name, condition, and value
that describe a named filter.
}

\usage{
## Constructors

\S4method{VcfFixedFilter}{character,character,ANY}(name, condition, value)

\S4method{VcfInfoFilter}{character,character,ANY}(name, condition, value)

\S4method{VcfVepFilter}{character,character,ANY}(name, condition, value)
}

\section{Accessors}{
In the following code snippets \code{x} is a \code{VcfBasicFilter} object.
\describe{
    \item{\code{name(x)}, \code{name(x) <- value}}{
    Gets or sets the name of the field to filter.
    }
    \item{\code{condition(x)}, \code{condition(x) <- value}}{
    Gets or sets the condition that values must pass, relative to \code{value}.
    }
    \item{\code{value(x)}, \code{value(x) <- value}}{
    Gets or sets the target values of the filter.
    }
    \item{\code{filterType(x)}}{
    Gets the type of filter.
    }
}
}

\arguments{
    \item{name}{
A \code{character} value of length 1 that defines the field to use for
filtering.
}
    \item{condition}{
A \code{character} value of length 1 that defines the binary operator used
to compare data to the target values.
}
    \item{value}{
Target value(s).
}
}

\value{
A \code{VcfBasicFilter} object.
}

\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\linkS4class{VcfFilterList}},
\code{\linkS4class{VCF}},
and \code{\link{parseCSQToGRanges,VCF-method}}.
}

\examples{
# Constructors ----

VcfFixedFilter(name = "FILTER", condition = "==", value = "PASS")

VcfInfoFilter(name = "MAF", condition = "<=", value = 1e-6)

filt <- VcfVepFilter(
    name = "IMPACT", condition = "\%in\%", value = c("MODERATE", "HIGH"))

# Accessors ----

## Name, condition, value, and type stored in the VcfBasicFilter object
name(filt)
condition(filt)
value(filt)
filterType(filt)
}
