\name{VcfBasicRules-class}
\docType{class}

% Class
\alias{class:VcfFixedRules}
\alias{VcfFixedRules-class}
\alias{class:VcfInfoRules}
\alias{VcfInfoRules-class}
\alias{class:VcfVepRules}
\alias{VcfVepRules-class}

% Constructors/Methods:
\alias{VcfFixedRules}
\alias{VcfInfoRules}
\alias{VcfVepRules}

\alias{VcfFixedRules}
\alias{VcfInfoRules}
\alias{VcfVepRules}

\alias{initialize,VcfFixedRules-method}
\alias{initialize,VcfInfoRules-method}
\alias{initialize,VcfVepRules-method}

% Accessor methods:
\alias{type,FilterRules-method}
\alias{type,VcfFixedRules-method}
\alias{type,VcfInfoRules-method}
\alias{type,VcfVepRules-method}

\alias{vep,FilterRules-method}
\alias{vep,VcfFixedRules-method}
\alias{vep,VcfInfoRules-method}
\alias{vep,VcfVepRules-method}
\alias{vep<-,VcfVepRules,character-method}

% Subsetting methods:
% [[
\alias{[[,VcfFixedRules,ANY,ANY-method}
\alias{[[,VcfInfoRules,ANY,ANY-method}
\alias{[[,VcfVepRules,ANY,ANY-method}
% [
\alias{[,VcfFixedRules,ANY,ANY-method}
\alias{[,VcfInfoRules,ANY,ANY-method}
\alias{[,VcfVepRules,ANY,ANY-method}

% Replacement methods:
% [[<-
\alias{[[<-,VcfFixedRules,ANY,ANY-method}
\alias{[[<-,VcfInfoRules,ANY,ANY-method}
\alias{[[<-,VcfVepRules,ANY,ANY-method}
% [<-
\alias{[<-,VcfFixedRules,numeric,missing,VcfFixedRules-method}
\alias{[<-,VcfInfoRules,numeric,missing,VcfInfoRules-method}
\alias{[<-,VcfVepRules,numeric,missing,VcfVepRules-method}

% Combine methods:
\alias{c,VcfFixedRules-method}
\alias{c,VcfInfoRules-method}
\alias{c,VcfVepRules-method}

\alias{append,VcfFixedRules,FilterRules-method}
\alias{append,VcfInfoRules,FilterRules-method}
\alias{append,VcfVepRules,FilterRules-method}

% Evaluation:
\alias{eval,VcfFixedRules,VCF-method}
\alias{eval,VcfInfoRules,VCF-method}
\alias{eval,VcfVepRules,VCF-method}

\title{VCF filters class objects sub-types}

\description{
The \code{VcfFixedRules} and \code{VcfInfoRules} classes
store filters applicable to the \code{fixed} and \code{info} slots
of \code{VCF} objects, respectively.

The \code{VcfVepRules} stores filters applicable to Ensembl VEP predictions
stores in a given INFO key.
}

\details{
All arguments are first passed to \code{S4Vectors::FilterRules}
before re-typing the resulting as a \code{VcfFixedRules}, \code{VcfInfoRules},
or \code{VcfVepRules} class.
}

\section{Accessor methods}{
In the following code snippets \code{x} is an object from any of the classes
decribed in this help page, except when specified otherwise.
\describe{
    \item{\code{active(x)}, \code{active(x)<-}}{
    Gets or sets the active state of each filter rule in \code{x}.
    Inherited from \code{\linkS4class{FilterRules}}
    }
    \item{\code{vep(x)}, \code{vep(x)<- }}{
    Gets or sets the INFO key where the Ensembl VEP predictions
    to use for filtering are stored.
    Returns \code{NA_character_} for filters not applicable to VEP predictions.
    }
    \item{\code{type(x)}}{
    Returns \code{"filter"} (\code{linkS4class{FilterRules}}),
    \code{"fixed"} (\code{linkS4class{VcfFixedRules}}),
    \code{"info"} (\code{linkS4class{VcfInfoRules}}),
    or \code{"vep"} (\code{linkS4class{VcfVepRules}})
    as a \code{character} vector of \code{length(x)}.
    }
}
}

\section{Constructors}{
\describe{}{
\code{VcfFixedRules(exprs = list(), ..., active = TRUE)}

\code{VcfInfoRules(exprs = list(), ..., active = TRUE)}

\code{VcfVepRules(exprs = list(), ..., active = TRUE, vep = "CSQ")}

All methods construct an object of the corresponding class
with the rules given in the list \code{exprs} or in \code{...}.
The initial active state of the rules is given by \code{active},
which is recycled as necessary.

See the constructor of \code{FilterRules} for more details.
}
}

\section{Subsetting and Replacement}{
In the following code snippets \code{x} and \code{value}
are objects from any of the classes described in this help page.

\describe{
    \item{}{
    \code{x[i]}: Subsets the filter rules using the same interface as for
    \code{\linkS4class{List}}.
    }
    \item{}{
    \code{x[[i]]}: Extracts an expression or function via the same interface
    as for \code{\linkS4class{List}}.
    }
    \item{}{
    \code{x[i] <- value}: Replaces a filter rule by one of the
    \strong{same} class.
    The active state(s) and name(s) are transferred from \code{value} to
    \code{x}.
    }
    \item{}{
    \code{x[[i]] <- value}:
    The same interface as for \code{\linkS4class{List}}.
    The default active state for new rules is TRUE.
    }
}
}

\section{Combining}{
In the following code snippets \code{x}, \code{values}, and \code{...}
are objects from any of the classes described in this help page, or
\code{VcfFilterRules}.

\describe{
    \item{}{
    \code{append(x, values, after = length(x))}:
    Appends the values onto \code{x} at the index given by \code{after}.
    }
    \item{}{
    \code{c(x, ...,)}:
    Concatenates the filters objects in \code{...} onto the end of \code{x}.
    }
}

Note that combining rules of different types
(\emph{e.g.} \code{VcfFixedRules} and \code{VcfVepRules})
produces a \code{VcfFilterRules} object.
}

\section{Evaluating}{
As described in the \code{S4Vectors} documentation:
\describe{
\item{}{
    \code{eval(expr, envir, enclos)}:
    Evaluates a rule instance (passed as the \code{expr} argument)
    in their respective context of a \code{VCF} object
    (passed as the \code{envir} argument).
    \emph{i.e.}:
    \itemize{
    \item{\strong{\code{VcfFixedRules}}: \code{fixed(envir)}}
    \item{\strong{\code{VcfInfoRules}}: \code{info(envir)}}
    \item{\strong{\code{VcfVepRules}}:
        \code{mcols(parseCSQToGRanges(envir, ...))}}
    \item{\strong{\code{FilterRules}}: \code{envir}}
    }
}
\item{}{
    \code{evalSeparately(expr, envir, enclos)}:

        \code{subsetByFilter(x, filter)}

        \code{summary(object)}

        See \code{\link{eval,FilterRules,ANY-method}} for details.
}
}
}

\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\linkS4class{FilterRules}},
\code{\linkS4class{VcfFilterRules}},
and \code{\linkS4class{VCF}}.
}

\examples{
# Constructors ----

fixedRules <- VcfFixedRules(list(
    pass = expression(FILTER == "PASS"),
    qual = expression(QUAL > 20)
    ))
fixedRules

infoRules <- VcfInfoRules(list(
    common = expression(MAF > 0.01), # minor allele frequency
    alt = expression(ALT > 0) # count of alternative homozygotes
    ))
infoRules

vepRules <- VcfVepRules(list(
    missense = expression(Consequence \%in\% c("missense_variant")),
    CADD = expression(CADD_PHRED > 15)
    ))
vepRules

filterRules <- S4Vectors::FilterRules(list(
    PASS = function(x) fixed(x)$FILTER == "PASS",
    COMMON = function(x) info(x)$MAF > 0.05
    ))
filterRules

# Accessors ----

## get/set the active state directly
S4Vectors::active(infoRules)
S4Vectors::active(infoRules)["common"] <- FALSE

## See S4Vectors::FilterRules for more examples


# Example data ----

# VCF file
vcfFile <- system.file("extdata", "moderate.vcf", package = "TVTB")

# TVTB parameters
tparam <- TVTBparam(Genotypes("0|0", c("0|1", "1|0"), "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(vcfFile, param = tparam)
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)
vcf <- addOverallFrequencies(vcf)


# Applying filters to VCF objects ----

## Evaluate filters
S4Vectors::eval(fixedRules, vcf)
S4Vectors::eval(infoRules, vcf)
S4Vectors::eval(vepRules, vcf)
S4Vectors::eval(filterRules, vcf)

summary(S4Vectors::eval(vepRules, vcf))


## Evaluate filters separately
S4Vectors::evalSeparately(vepRules, vcf)

summary(S4Vectors::evalSeparately(vepRules, vcf))

## Subset VCF by filters
S4Vectors::subsetByFilter(vcf, vepRules)


# Subsetting and Replacement ----

vep1 <- vepRules[1] # VcfVepRules
vepRules[[1]] # expression

# Replace the expression (active reset to TRUE, original name retained)
vepRules[[2]] <- expression(CADD_PHRED > 30)

# Replace the rule (active state and name transferred from v5obj)
vepRules[2] <- VcfVepRules(
    list(newRule = expression(CADD_PHRED > 30)),
    active = FALSE)
}
