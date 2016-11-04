\name{TVTBparam-class}
\docType{class}

% Class
\alias{class:TVTBparam}
\alias{TVTBparam-class}

% Constructors/Methods:
\alias{TVTBparam}
\alias{TVTBparam-methods}
\alias{TVTBparam,list-method}
\alias{initialize,TVTBparam-method}

% Accessor methods:
\alias{genos}
\alias{genos,TVTBparam-method}
\alias{genos<-}
\alias{genos<-,TVTBparam,list-method}
\alias{ranges}
\alias{ranges,TVTBparam-method}
\alias{ranges<-}
\alias{ranges<-,TVTBparam,GRangesList-method}
\alias{aaf}
\alias{aaf,TVTBparam-method}
\alias{aaf<-}
\alias{aaf<-,TVTBparam,character-method}
\alias{maf}
\alias{maf,TVTBparam-method}
\alias{maf<-}
\alias{maf<-,TVTBparam,character-method}
\alias{vep}
\alias{vep,TVTBparam-method}
\alias{vep<-}
\alias{vep<-,TVTBparam,character-method}
\alias{bp}
\alias{bp,TVTBparam-method}
\alias{bp<-}
\alias{bp<-,TVTBparam,BiocParallelParam-method}

% Other methods / functions
\alias{hRef}
\alias{hRef,TVTBparam-method}
\alias{hRef<-}
\alias{hRef<-,TVTBparam,list-method}
\alias{hRef<-,TVTBparam,character-method}
\alias{het}
\alias{het,TVTBparam-method}
\alias{het<-}
\alias{het<-,TVTBparam,list-method}
\alias{het<-,TVTBparam,character-method}
\alias{hAlt}
\alias{hAlt,TVTBparam-method}
\alias{hAlt<-}
\alias{hAlt<-,TVTBparam,list-method}
\alias{hAlt<-,TVTBparam,character-method}
\alias{carrier}
\alias{carrier,TVTBparam-method}
\alias{carrier<-}
\alias{carrier<-,TVTBparam,list-method}

\title{TVTBparam class objects}

\description{
The \code{TVTBparam} class stores recurrent parameters of the \code{TVTB}
package in a convenient format.
}

\usage{
\S4method{TVTBparam}{list}(
    genos,
    ranges = GRangesList(),
    aaf = "AAF", maf = "MAF", vep = "CSQ", bp = SerialParam())
}

\section{Accessor methods}{
In the following code snippets \code{x} is a \code{TVTBparam} object.
\describe{
    \item{\code{genos(x)}, \code{genos(x) <- value}}{
    Gets or sets a named \code{list} of length three
    containing \code{character} vectors that declare
    homozygote reference, heterozygote, and homozygote alternate genotypes,
    respectively.
    If a named \code{list}, the names declare the individual suffixes
    used to store the counts of each set of genotypes in the INFO field of the
    \code{VCF} object.
    See also \code{hRef}, \code{het}, \code{hAlt}, and \code{carrier}
    accessors.
    }
    \item{\code{ranges(x)}, \code{ranges(x) <- value}}{
    Gets or sets a \code{GRangesList}
    that declares a list of genomic ranges to import and group variants.
    }
    \item{\code{hRef(x)}, \code{hRef(x) <- value}}{
    Gets or sets the (named) \code{list} of length 1
    or unnamed \code{character} vector
    that declare homozygote reference genotypes.
    If a named \code{list}, the name declares the suffix used to store the
    counts of the set of homozygote reference genotypes in the INFO field of
    the \code{VCF} object.
    }
    \item{\code{het(x)}, \code{het(x) <- value}}{
    Gets or sets the (named) \code{list} of length 1
    or unnamed \code{character} vector
    that declare heterozygote genotypes.
    If a named \code{list}, the name declares the suffix used to store the
    counts of the set of homozygote heterozygote genotypes in the INFO field
    of the \code{VCF} object.
    }
    \item{\code{hAlt(x)}, \code{hAlt(x) <- value}}{
    Gets or sets the (named) \code{list} of length 1
    or unnamed \code{character} vector
    that declare homozygote alternate genotypes.
    If a named \code{list}, the name declares the suffix used to store the
    counts of the set of homozygote alternate genotypes in the INFO field of
    the \code{VCF} object.
    }
    \item{\code{carrier(x)}, \code{carrier(x) <- value}}{
    Gets or sets a named \code{list} of length two
    containing \code{character} vectors that declare
    heterozygote, and homozygote alternate genotypes,
    respectively.
    If a named \code{list}, the names declare the individual suffixes
    used to store the counts of each set of genotypes in the INFO field of the
    \code{VCF} object.
    }
    \item{\code{aaf(x)}, \code{aaf(x) <- value}}{
    Gets or sets the INFO key suffix used to store the
    alternate allele frequency (AAF).
    }
    \item{\code{maf(x)}, \code{maf(x) <- value}}{
    Gets or sets the INFO key suffix used to store the
    minor allele frequency (MAF).
    }
    \item{\code{vep(x)}, \code{maf(x) <- value}}{
    Gets or sets the INFO key suffix used to extract the
    VEP predictions.
    }
    \item{\code{bp(x)}, \code{bp(x) <- value}}{
    Gets or sets the \code{BiocParallel} parameters.
    }
}
}

\arguments{
    \item{genos}{
A \code{list} of length three
containing \code{character} vectors that declares the three sets of
homozygote reference, heterozygote, and homozygote alternate genotypes,
in this order.

If named, \code{names} set the individual INFO key suffix used to
store the counts of each set of genotypes.
If unnamed, \code{c("REF","HET","ALT")} are used as default key
suffixes.
See also \emph{Details} section.
}
    \item{ranges}{
A \code{GRangesList} of genomic regions.
\emph{In the future, may be used to facet statistics and figures.}
}
\item{aaf}{
INFO key suffix used to store the alternate allele frequency (AAF).
}
    \item{maf}{
INFO key suffix used to store the minor allele frequency (MAF).
}
    \item{vep}{
INFO key suffix used to extract the VEP predictions.
}
    \item{bp}{A \code{BiocParallelParam} object.
}
}

\value{
A \code{TVTBparam} object that contains recurrent parameters.
}

\details{
For each \emph{suffix} stored in the \code{TVTBparam} object, \code{TVTB}
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
\code{\linkS4class{ExpandedVCF}},
\code{\link{addCountGenos-methods}}
\code{\link{vepInPhenoLevel-methods}},
\code{\link{densityVep-methods}},
\code{\link{tabulateVep-methods}},
\code{\link{variantsInSamples-methods}},
and \code{\linkS4class{BiocParallelParam}}.
}

\examples{
# Constructors ----

genotypes <- list(
    REF = c("0|0"),
    HET = c("0|1", "1|0"),
    ALT = c("1|1")
)

grl <- GenomicRanges::GRangesList(GenomicRanges::GRanges(
    seqnames = "15",
    ranges = IRanges::IRanges(
        start = 48413170, end = 48434757, names = "SLC24A5")
    ))

tparam <- TVTBparam(genos = genotypes, ranges = grl)

# Accessors ----

## Homozygote reference, heterozygote, and alternate heterozygote genotypes
## stored in the TVTBparam object returned by the genos() accessor.
genos(tparam)

## Genomic ranges stored in the TVTBparam object returned by the ranges()
## accessor.
ranges(tparam)

## Individual genotypes can be extracted with hRef(), het(), hAlt() accessors.
hRef(tparam)
het(tparam)
hAlt(tparam)

## Their individual INFO key suffixes can be extracted with names() applied to
## the above accessors.
names(hRef(tparam))
names(het(tparam))
names(hAlt(tparam))

## Heterozygote, and alternate heterozygote genotypes are
## returned by the carrier() accessor.
carrier(tparam)
names(carrier(tparam))

## INFO key suffix of alternate/minor allele frequency returned by the aaf()
## and maf() accessors.
aaf(tparam)
maf(tparam)

## INFO key suffix of the VEP predictions returned by the vep() accessor.
vep(tparam)

## BiocParallel parameters
bp(tparam)
}
