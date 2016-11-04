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
\alias{vep,VcfVepRules-method}

% Subsetter methods:
\alias{[[<-,VcfFixedRules,ANY,ANY-method}
\alias{[[<-,VcfInfoRules,ANY,ANY-method}
\alias{[[<-,VcfVepRules,ANY,ANY-method}

% Combine methods:
\alias{c,VcfFixedRules-method}
\alias{c,VcfInfoRules-method}
\alias{c,VcfVepRules-method}

% Evaluation:
\alias{eval,VcfFixedRules,ExpandedVCF-method}
\alias{eval,VcfInfoRules,ExpandedVCF-method}
\alias{eval,VcfVepRules,ExpandedVCF-method}
\alias{eval,VcfFilterRules,ExpandedVCF-method}

\title{VCF filters class objects sub-types}

\description{
The \code{VcfFixedRules} and \code{VcfInfoRules} classes
store filters applicable to the \code{fixed} and \code{info} slots
of the \code{ExpandedVCF} objects, respectively.

The \code{VcfVepRules} stores filters applicable to Ensembl VEP predictions
stores in a given INFO key.
}

\details{
All arguments are first passed to \code{S4Vectors::FilterRules}
before re-typing the resulting as a \code{VcfFixedRules}, \code{VcfInfoRules},
or \code{VcfVepRules} class.
}

\section{Accessor methods}{
In the following code snippets \code{x} is a \code{VcfFilterRules} object.
\describe{
    \item{\code{vep(x)}}{
    Gets the INFO key where the Ensembl VEP predictions to use for filtering
    are stored.
    }
}
And accessors inherited from \code{S4Vectors::FilterRules} class.
}

\section{Constructors}{
\describe{}{
\code{VcfFixedRules(exprs = list(), ..., active = TRUE)}

\code{VcfInfoRules(exprs = list(), ..., active = TRUE)}

\code{VcfVepRules(exprs = list(), ..., active = TRUE, vep = "CSQ")}

all construct an object of the corresponding class
with the rules given in the list \code{exprs} or in \code{...}.
The initial active state of the rules is given by \code{active},
which is recycled as necessary.

See the constructor of \code{FilterRules} for more details.
}
}

\section{Subsetting and Replacement}{
Inherited from \code{S4Vectors::FilterRules} class.
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
    Evaluates a rule instance (passed as the \code{expr} argument)
    in their respective context of a \code{ExpandedVCF} object
    (passed as the \code{envir} argument).
    \emph{I.e.}:
    \itemize{
    \item{\strong{\code{VcfFixedRules}}: \code{fixed(envir)}}
    \item{\strong{\code{VcfInfoRules}}: \code{info(envir)}}
    \item{\strong{\code{VcfVepRules}}: \code{parseCSQToGRanges(envir, ...)}}
    }

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
    rule instance (passed as the \code{expr}
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
\code{\linkS4class{FilterRules}},
\code{\linkS4class{VcfFilterRules}},
and \code{\linkS4class{ExpandedVCF}}.
}

\examples{
# Constructors ----

fixedRules <- VcfFixedRules(exprs = list(
    pass = expression(FILTER == "PASS"),
    qual = expression(QUAL > 20)
    ))
fixedRules

infoRules <- VcfInfoRules(exprs = list(
    common = expression(MAF > 0.1), # minor allele frequency
    alt = expression(ALT > 0) # count of alternative homozygotes
    ))
infoRules

vepRules <- VcfVepRules(exprs = list(
    missense = expression(Consequence \%in\% c("missense_variant")),
    CADD = expression(CADD_PHRED > 15)
    ))
vepRules

# Accessors ----

## get/set the active state directly
S4Vectors::active(infoRules)
S4Vectors::active(infoRules)["common"] <- FALSE

## See S4Vectors::FilterRules for more examples

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
vcf <- VariantAnnotation::expand(vcf, row.names = TRUE)
vcf <- addOverallFrequencies(vcf = vcf, param = tparam)

# Applying filters to ExpandedVCF objects ----

## Evaluate filters
S4Vectors::eval(expr = fixedRules, envir = vcf)
S4Vectors::eval(expr = infoRules, envir = vcf)
S4Vectors::eval(expr = vepRules, envir = vcf)
summary(S4Vectors::eval(expr = vepRules, envir = vcf))

## Evaluate filters separately
S4Vectors::evalSeparately(expr = vepRules, envir = vcf)
summary(S4Vectors::evalSeparately(expr = vepRules, envir = vcf))

## Subset ExpandedVCF by filters
S4Vectors::subsetByFilter(x = vcf, filter = vepRules)

# Subsetting and Replacement ----

vep1 <- vepRules[1] # VcfVepRules
vepRules[[1]] # expression

vepRules[[2]] <- expression(CADD_PHRED > 30)
}
