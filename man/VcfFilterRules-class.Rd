\name{VcfFilterRules-class}
\docType{class}

% Class
\alias{class:VcfFilterRules}
\alias{VcfFilterRules-class}

% Constructors/Methods:
\alias{VcfFilterRules}
%\alias{VcfFilterRules,ANY-method}
\alias{initialize,VcfFilterRules-method}

% Accessor methods:
\alias{type}
\alias{type,VcfFilterRules-method}
\alias{vep,VcfFilterRules-method}

% Subsetter methods:
\alias{[,VcfFilterRules,ANY,ANY,missing-method}
\alias{[,VcfFilterRules,ANY,ANY,logical-method}
\alias{[[,VcfFilterRules,ANY,ANY-method}
\alias{[[<-,VcfFilterRules,ANY,ANY-method}

% Combine methods:
\alias{c,VcfFilterRules-method}

% Evaluation:
\alias{eval-methods}

\title{VcfFilterRules class objects}

\description{
The \code{VcfFilterRules} class stores various types of filters applicable to
\code{ExpandedVCF} objects.
}

\details{
All arguments must be \code{VcfFixedRules}, \code{VcfInfoRules},
or \code{VcfVepRules} objects.
}

\section{Accessor methods}{
In the following code snippets \code{x} is a \code{VcfFilterRules} object.
\describe{
    \item{\code{type(x)}}{
    Gets the type of each filter stored in a \code{VcfFilterRules} object.
    }
    \item{\code{vep(x)}}{
    Gets the INFO key where the Ensembl VEP predictions to use for filtering
    are stored.
    }
}
And accessors inherited from \code{S4Vectors::FilterRules} class.
}

\section{Constructors}{
\describe{}{
\code{VcfFilterRules(...)}
constructs an \code{VcfFilterRules} object from
\code{VcfFixedRules}, \code{VcfInfoRules}, and \code{VcfVepRules} objects
in \code{...}.
}
}

\section{Subsetting and Replacement}{
In the code snippets below, \code{x} is a \code{VcfFilterRules} object.

\describe{
\item{}{
    \code{x[i]}: Subsets the filter rules using the
    same interface as for \code{\linkS4class{Vector}}.
    If all filter rules are of the same type, re-type the resulting
    filter rules accordingly.
}
\item{}{
    \code{x[[i]]}: Extracts a
    \code{VcfFixedRules}, \code{VcfInfoRules}, or \code{VcfVepRules} object
    via the same interface as for \code{\linkS4class{List}}.
}
\item{}{
    \code{x[[i]] <- value}: The same interface as for
    \code{\linkS4class{List}}. The default active state for new
    rules is \code{TRUE}.
}
}
}

\section{Combining}{
Inherited from \code{S4Vectors::FilterRules} class.
}

\section{Evaluating}{
\describe{
\item{}{
    \code{eval(expr, envir = parent.frame(),
    enclos = if (is.list(envir) || is.pairlist(envir))
    parent.frame() else baseenv())}:
    Evaluates each active rule in a \code{VcfFilterRules} instance
    (passed as the \code{expr} argument)
    in their respective context of a \code{ExpandedVCF} object
    (passed as the \code{envir} argument).

    Expression rules are
    evaluated in \code{envir}, while function rules are invoked with
    \code{envir} as their only argument. The evaluation of a rule
    should yield a logical vector. The results from the rule
    evaluations are combined via the AND operation (i.e. \code{&}) so
    that a single logical vector is returned from \code{eval}.
}
\item{}{
    \code{evalSeparately(expr, envir = parent.frame(), enclos = if
    (is.list(envir) || is.pairlist(envir)) parent.frame() else
    baseenv())}: Evaluates separately each rule in a
    \code{VcfFilterRules} instance (passed as the \code{expr}
    argument). Expression rules are evaluated in \code{envir}, while
    function rules are invoked with \code{envir} as their only
    argument. The evaluation of a rule should yield a logical
    vector. The results from the rule evaluations are combined into
    a logical matrix, with a column for each rule. This is
    essentially the parallel evaluator, while \code{eval} is the
    serial evaluator.
}
\item{}{
    \code{subsetByFilter(x, filter)}: Evaluates \code{filter} on
    \code{x} and uses the result to subset \code{x}. The result
    contains only the elements in \code{x} for which \code{filter}
    evaluates to \code{TRUE}.
}
\item{}{\code{summary(object)}:
    Returns an integer vector with the number of elements
    that pass each rule in a \code{FilterMatrix}
    (passed as the \code{object} argument), along with
    a count of the elements that pass all filters.
}
}
}

\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\linkS4class{VcfFixedRules}},
\code{\linkS4class{VcfInfoRules}},
\code{\linkS4class{VcfVepRules}},
and \code{\linkS4class{ExpandedVCF}}.
}

\examples{
# Constructors ----

fixedR <- VcfFixedRules(exprs = list(
    pass = expression(FILTER == "PASS"),
    qual = expression(QUAL > 20)
    ))
fixedR

infoR <- VcfInfoRules(exprs = list(
    common = expression(MAF > 0.1), # minor allele frequency
    present = expression(ALT + HET > 0) # count of non-REF homozygotes
    ))
# ...is synonym to...
infoR <- VcfInfoRules(exprs = list(
    common = expression(MAF > 0.1), # minor allele frequency
    present = expression(ALT > 0 | HET > 0)
    ))
infoR

vepR <- VcfVepRules(exprs = list(
    missense = expression(Consequence \%in\% c("missense_variant")),
    CADD = expression(CADD_PHRED > 15)
    ))
vepR

vcfRules <- VcfFilterRules(fixedR, infoR, vepR)
vcfRules

# Accessors ----

## Type of each filter stored in the VcfFilterRules object
Biostrings::type(vcfRules)

# Example data ----

# VCF file
extdata <- file.path(system.file(package = "TVTB"), "extdata")
vcfFile <- file.path(extdata, "moderate.vcf")

# TVTB parameters
tparam <- TVTBparam(
    genos = list(
        REF = "0|0",
        HET = c("0|1", "1|0"),
        ALT = "1|1"))

# Pre-process variants
vcf <- VariantAnnotation::readVcf(file = vcfFile)
vcf <- VariantAnnotation::expand(vcf)
vcf <- addOverallFrequencies(vcf = vcf, param = tparam)

# Applying filters to ExpandedVCF objects ----

## Evaluate filters
S4Vectors::eval(expr = vcfRules, envir = vcf)

## Evaluate filters separately
as.data.frame(S4Vectors::evalSeparately(expr = vcfRules, envir = vcf))

# Interestingly, the only common missense variant has a lower CADD score
## Deactivate the CADD score filter
S4Vectors::active(vcfRules)["CADD"] <- FALSE

## Subset ExpandedVCF by filters (except CADD, deactivated above)
S4Vectors::subsetByFilter(x = vcf, filter = vcfRules)

# Subsetting and Replacement ----

v123 <- vcfRules[1:3]
v5 <- vcfRules[[5]]
v123[[2]] <- v5
}
