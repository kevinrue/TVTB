\name{NamedFilter-class}
\docType{class}

% Class
\alias{class:NamedFilter}
\alias{NamedFilter-class}

% Constructors/Methods:
\alias{NamedFilter}
\alias{NamedFilter-methods}
\alias{NamedFilter,character,character-method}
\alias{NamedFilter,character,character,ANY-method}
\alias{initialize,NamedFilter-method}

% Accessor methods:
\alias{name}
\alias{name,NamedFilter-method}
\alias{name<-}
\alias{name<-,NamedFilter,character-method}
\alias{condition}
\alias{condition,NamedFilter-method}
\alias{condition<-}
\alias{condition<-,NamedFilter,character-method}
\alias{value}
\alias{value,NamedFilter-method}
\alias{value<-}
\alias{value<-,NamedFilter-method}

% Other methods / functions

\title{NamedFilter class objects}

\description{
The \code{NamedFilter} class stores the name, condition, and value
that describe a named filter.
}

\usage{
## Constructors

\S4method{NamedFilter}{character,character,ANY}(name, condition, value)
}

\section{Accessors}{
In the following code snippets \code{x} is a \code{NamedFilter} object.
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
}
}

\arguments{
    \item{name}{
A \code{character} value of length 1 that defines the field to filter.
}
    \item{condition}{
A \code{character} value of length 1 that defines the binary operator used
to compare data to the target values.
}
    \item{value}{
Target values.
}
}

\value{
A \code{NamedFilter} object that contains recurrent parameters.
}

\author{
Kevin Rue-Albrecht
}

%\seealso{
% Soon: NamedFilterList, VcfFilter
%}

\examples{
# Constructors ----

NamedFilter(name = "FILTER", condition = "==", value = "PASS")

NamedFilter(name = "MAF", condition = "<=", value = 1e-6)

nf <- NamedFilter(
    name = "IMPACT", condition = "\%in\%", value = c("MODERATE", "HIGH"))

# Accessors ----

## Name, condition, and value stored in the tSVEParam object
name(nf)
condition(nf)
value(nf)
}
