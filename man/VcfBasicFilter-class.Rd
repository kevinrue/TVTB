\name{VcfBasicFilter-class}
\docType{class}

% Class
\alias{class:VcfBasicFilter}
\alias{VcfBasicFilter-class}
\alias{class:VcfInfoFilter}
\alias{VcfInfoFilter-class}

% Constructors/Methods:
\alias{VcfInfoFilter}
\alias{VcfInfoFilter-methods}
\alias{VcfInfoFilter,character,character-method}
\alias{VcfInfoFilter,character,character,ANY-method}
\alias{initialize,VcfInfoFilter-method}

% Accessor methods:
\alias{name}
\alias{name,VcfInfoFilter-method}
\alias{name<-}
\alias{name<-,VcfInfoFilter,character-method}
\alias{condition}
\alias{condition,VcfInfoFilter-method}
\alias{condition<-}
\alias{condition<-,VcfInfoFilter,character-method}
\alias{value}
\alias{value,VcfInfoFilter-method}
\alias{value<-}
\alias{value<-,VcfInfoFilter-method}

% Other methods / functions

\title{VcfBasicFilter class objects}

\description{
The \code{VcfInfoFilter} class stores the name, condition, and value
that describe a named filter.
}

\usage{
## Constructors

\S4method{VcfInfoFilter}{character,character,ANY}(name, condition, value)
}

\section{Accessors}{
In the following code snippets \code{x} is a \code{VcfInfoFilter} object.
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
A \code{VcfInfoFilter} object that contains recurrent parameters.
}

\author{
Kevin Rue-Albrecht
}

%\seealso{
% Soon: VcfInfoFilterList, VcfFilter
%}

\examples{
# Constructors ----

VcfInfoFilter(name = "FILTER", condition = "==", value = "PASS")

VcfInfoFilter(name = "MAF", condition = "<=", value = 1e-6)

nf <- VcfInfoFilter(
    name = "IMPACT", condition = "\%in\%", value = c("MODERATE", "HIGH"))

# Accessors ----

## Name, condition, and value stored in the tSVEParam object
name(nf)
condition(nf)
value(nf)
}
