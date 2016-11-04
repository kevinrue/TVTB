\name{TVTB-NEWS}
\title{NEWS: The VCF Tool Box}

\section{Changes in version 0.99.1 2016-10-10}{
\subsection{Major changes}{
\itemize{
\item
    New dedicated \code{Genotypes} class to store homozygote reference,
    heterozygote, and homozygote alternate genotype codes,
    along with the suffixes that define the INFO keys
    used to store their respective data in the VCF object.
\item
    Removed families of methods \code{tabulate*} and \code{density*} from the
    \code{NAMESPACE}. The features may be revisited in the future.
    The associated code and documentation was saved in the \code{inst/sandbox}
    subfolder for future reference.
\item
    New slot \code{svp} in \code{TVTBparam} class to store \code{ScanVcfParam}
    objects. New associated accessor methods.
    Moreover, \code{TVTBparam} may be coerced to \code{ScanVcfParam}.
\item
    New signatures for method \code{readVcf} that supports
    \code{param=TVTBparam}, and optional phenotypes.
    The method stores \code{TVTBparam} in the \code{metadata} slot of
    the \code{VCF} object, and phenotypes using the \code{colData} accessor.
\item
    \code{TVTBparam} are no longer an argument of downstream methods; instead,
    they must be stored in \code{metadata(vcf)[["TVTBparam"]]}
}
}
\subsection{Minor changes}{
\itemize{
\item
    \code{hRef} and \code{hAlt} accessors renamed to \code{ref}and \code{alt},
    respectively.
\item
    \code{suffix} accessor to returned named \code{character} vector
    for classes \code{Genotypes} and \code{TVTBparam}.
\item
    Removed functions relevant only to the \emph{Shiny} application from the
    \code{NAMESPACE} (\code{getEdb}, \code{EnsDbFilter}, \code{chr2file}).
\item
    \code{VcfFilterRules} can also store instances of the parent
    \code{FilterRules} class.
\item
    Defined default return value for accessors \code{vep} and \code{type},
    to avoid unnecessary \code{switch} statements.
\item
    Simplified code of inherited methods following updates to the relevant
    packages (\emph{e.g.} \code{S4Vectors}).
\item
    Better respect of coding standards: removed superfluous usage of
    explicit argument naming, removed superfluous \code{initialize} methods.
\item
    Man pages, vignettes, unit tests and \emph{Shiny} application updated
    to reflect changes to the package.
}
}
}

\section{Changes in version 0.99.0 2015-09-15}{
\subsection{New features}{
\itemize{
\item
    First release submitted to the Bioconductor review process
    See \code{DESCRIPTION} file for details.
}
}
}
