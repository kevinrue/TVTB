\name{Genotypes-class}
\docType{class}

% Class
\alias{class:Genotypes}
\alias{Genotypes-class}

% Constructors/Methods:
\alias{Genotypes}
\alias{Genotypes-methods}
\alias{Genotypes,character,character,character-method}
\alias{initialize,Genotypes-method}

% Accessor methods:
\alias{genos}
\alias{genos,Genotypes-method}

% Other methods / functions
\alias{ref,Genotypes-method}
\alias{ref<-,Genotypes,character-method}
\alias{het}
\alias{het,Genotypes-method}
\alias{het<-}
\alias{het<-,Genotypes,character-method}
\alias{alt,Genotypes-method}
\alias{alt<-,Genotypes,character-method}
\alias{carrier}
\alias{carrier,Genotypes-method}
\alias{suffix}
\alias{suffix,Genotypes-method}

\title{Genotypes class objects}

\description{
The \code{Genotypes} class stores genotype definitions in a convenient format.
}

\usage{
\S4method{Genotypes}{character,character,character}(
    ref, het, alt, suffix = c(ref="REF", het="HET", alt="ALT"))
}

\section{Accessor methods}{
In the following code snippets \code{x} is a \code{Genotypes} object.
\describe{
    \item{\code{ref(x)}, \code{ref(x) <- value}}{
    Gets or sets the \code{character} vector
    that declares homozygote reference genotypes.
    }
    \item{\code{het(x)}, \code{het(x) <- value}}{
    Gets or sets the \code{character} vector
    that declares heterozygote genotypes.
    }
    \item{\code{alt(x)}, \code{alt(x) <- value}}{
    Gets or sets the \code{character} vector
    that declares homozygote alternate genotypes.
    }
    \item{\code{genos(x)}}{
    Gets a \code{character} vector of concatenated
    homozygote reference, heterozygote, and homozygote alternate genotypes.
    See also \code{ref}, \code{het}, \code{alt}, and \code{carrier}
    accessors.
    }
    \item{\code{carrier(x)}}{
    Gets a \code{character} vectors of concatenated
    heterozygote and homozygote alternate genotypes.
    See also \code{het} and \code{alt} accessors.
    }
    \item{\code{suffix(x)}}{
    Gets a named \code{character} vector that declares individual suffixes
    used to store the data for each set of genotypes in the INFO field of the
    \code{VCF} object.
    Names of this vector are \code{ref}, \code{het}, and \code{alt}.
    }
}
}

\arguments{
    \item{ref}{
A \code{character} vector declaring the encoding of
homozygote reference genotypes.
}
    \item{het}{
A \code{character} vector declaring the encoding of
heterozygote genotypes.
}
    \item{alt}{
A \code{character} vector declaring the encoding of
homozygote alternate genotypes.
}
    \item{suffix}{
Set the individual INFO key suffixes used to
store the statistics of
homozygote reference, heterozygote, and homozygote alternate genotypes,
in this order.
See \emph{Details} section.
}
}

\value{
A \code{Genotypes} object that contains genotype definitions.
}

\details{
For each \emph{suffix} stored in the \code{Genotypes} object, \code{TVTB}
may store data in the \code{VCF} object under the INFO keys defined as follows:
\describe{
    \item{\emph{suffix}}{
    Statistics across all samples in the \code{ExpandedVCF}
    (\emph{e.g.} "MAF").}
    \item{\emph{phenotype_level_suffix}}{
    Statistics across samples associated with a given level of a given
    phenotype (\emph{e.g.} "gender_male_MAF").}
    }

Users are recommended to avoid using those INFO keys for other purposes.
}

\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\linkS4class{VCF}},
\code{\linkS4class{TVTBparam}},
and \code{\link{addCountGenos-methods}}.
}

\examples{
# Constructors ----

genotypes <- Genotypes("0|0", c("0|1", "1|0"), "1|1")


# Accessors ----

## Concatenated homozygote reference, heterozygote, and alternate heterozygote
## genotypes stored in the Genotypes object returned by the genos() accessor.
genos(genotypes)

## Individual genotypes can be extracted with ref(), het(), alt() accessors.
ref(genotypes)
het(genotypes)
alt(genotypes)

## Their individual INFO key suffixes can be extracted with suffix() accessors
## and the relevant name
suffix(genotypes)
suffix(genotypes)["ref"]
suffix(genotypes)["het"]
suffix(genotypes)["alt"]

## Concatenated heterozygote, and alternate heterozygote genotypes are
## returned by the carrier() accessor.
carrier(genotypes)
names(carrier(genotypes))
}
