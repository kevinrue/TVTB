\name{TVTBparam-class}
\docType{class}

% Class
\alias{class:TVTBparam}
\alias{TVTBparam-class}

% Constructors/Methods:
\alias{TVTBparam}
\alias{TVTBparam-methods}
\alias{TVTBparam,Genotypes-method}
\alias{initialize,TVTBparam-method}

% Accessor methods:
\alias{genos,TVTBparam-method}
\alias{genos<-}
\alias{genos<-,TVTBparam,Genotypes-method}
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
\alias{svp}
\alias{svp,TVTBparam-method}
\alias{svp<-}
\alias{svp<-,TVTBparam,ScanVcfParam-method}

% Other methods / functions
\alias{ref,TVTBparam-method}
\alias{ref<-,TVTBparam,list-method}
\alias{ref<-,TVTBparam,character-method}
\alias{het,TVTBparam-method}
\alias{het<-,TVTBparam,list-method}
\alias{het<-,TVTBparam,character-method}
\alias{alt,TVTBparam-method}
\alias{alt<-,TVTBparam,list-method}
\alias{alt<-,TVTBparam,character-method}
\alias{carrier,TVTBparam-method}
\alias{suffix,TVTBparam-method}

\title{TVTBparam class objects}

\description{
The \code{TVTBparam} class stores recurrent parameters of the \code{TVTB}
package in a convenient format.
}

\usage{
\S4method{TVTBparam}{Genotypes}(
        genos, ranges = GRangesList(),
        aaf = "AAF", maf = "MAF", vep = "CSQ", bp = SerialParam(),
        svp = ScanVcfParam(which = reduce(unlist(ranges))))
}

\section{Accessor methods}{
In the following code snippets \code{x} is a \code{TVTBparam} object.
\describe{
    \item{\code{genos(x)}, \code{genos(x) <- value}}{
    Gets or sets the \code{Genotypes} object stored in the \code{genos} slot.
    }
    \item{\code{ranges(x)}, \code{ranges(x) <- value}}{
    List of genomic ranges to group variants during analyses and plots.
    }
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
    \item{\code{carrier(x)}}{
    Gets a \code{character} vectors of concatenated
    heterozygote and homozygote alternate genotypes.
    See also \code{het} and \code{alt} accessors.
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
    \item{\code{suffix(x)}}{
    Gets a named \code{character} vector that declares individual suffixes
    used to store the data for each set of genotypes in the INFO field of the
    \code{VCF} object.
    Names of this vector are \code{ref}, \code{het}, \code{alt}, \code{aaf},
    and \code{maf}.
    }
    \item{\code{svp(x)}, \code{svp(x) <- value}}{
    Gets or sets the \code{ScanVcfParam} parameters.
    }
}
}

\arguments{
    \item{genos}{
A \code{Genotypes} object that declares the three sets of
homozygote reference, heterozygote, and homozygote alternate genotypes,
as well as the individual key suffix used to
store data for each set of genotypes in the \code{info} slot of a
\code{VCF} object.
See also \emph{Details} section.
}
    \item{ranges}{
A \code{GRangesList} of genomic regions. See \code{svp} argument.
\emph{In the future, may be used to facet statistics and figures.}
}
\item{aaf}{
INFO key suffix used to store the alternate allele frequency (AAF).
}
    \item{maf}{
INFO key suffix used to store the minor allele frequency (MAF).
}
    \item{vep}{
INFO key suffix used to extract the VEP predictions. See \code{svp} argument.
}
    \item{bp}{A \code{BiocParallelParam} object.
}
    \item{svp}{A \code{ScanVcfParam} object.
If none is supplied, the \code{\linkS4class{ScanVcfParam}} slot \code{which}
is automatically set to \code{reduce(unlist(ranges))}.
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
\code{\linkS4class{Genotypes}},
\code{\linkS4class{VCF}},
\code{\linkS4class{ExpandedVCF}},
\code{\link{addCountGenos-methods}}
\code{\link{vepInPhenoLevel-methods}},
\code{\link{variantsInSamples-methods}},
and \code{\linkS4class{BiocParallelParam}}.
}

\examples{
# Constructors ----

grl <- GenomicRanges::GRangesList(GenomicRanges::GRanges(
    "15", IRanges::IRanges(48413170, 48434757, names = "SLC24A5")
    ))

tparam <- TVTBparam(Genotypes("0|0", c("0|1", "1|0"), "1|1"), ranges = grl)


# Accessors ----

## The Genotypes object stored in the genos slot of the TVTBparam object
## return by the genos() accessor.
genos(tparam)

## Genomic ranges stored in the TVTBparam object returned by the ranges()
## accessor.
ranges(tparam)

## Individual genotypes can be extracted with ref(), het(), alt() accessors.
ref(tparam)
het(tparam)
alt(tparam)

## Their individual INFO key suffixes can be extracted with suffix() applied to
## the above accessors.
suffix(tparam)
suffix(tparam)["ref"]
suffix(tparam)["het"]
suffix(tparam)["alt"]
suffix(tparam)["aaf"]
suffix(tparam)["maf"]

## Heterozygote, and alternate heterozygote genotypes are
## returned by the carrier() accessor.
carrier(tparam)

## INFO key suffix of alternate/minor allele frequency returned by the aaf()
## and maf() accessors.
aaf(tparam)
maf(tparam)

## INFO key suffix of the VEP predictions returned by the vep() accessor.
vep(tparam)

## BiocParallel parameters
bp(tparam)

## ScanVcfParam parameters
svp(tparam)
}
