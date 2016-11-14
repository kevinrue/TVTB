\name{TVTB-NEWS}
\title{NEWS: The VCF Tool Box}

\section{Changes in version 1.1.8 2016-11-14}{
\subsection{New features}{
\itemize{
\item New method \code{pairsInfo} to visualise a matrix of pairwise plots
    that displays a metric calculated in levels of a given phenotype,
    and stored in columns of the \code{info} slot of a \code{VCF} object.
\item \code{ggpairs} method imported from the \code{GGally} package.
}
}
\subsection{Bug fix}{
\itemize{
\item Updated pattern used to detect INFO columns for a given metric, to be
    more specific.
}
}
\subsection{Major changes}{
\itemize{
\item Internal method \code{.findInfoMetricColumns} moved to \code{utils.R},
    as it is now used in two different user-visible methods.
\item Updated \code{Introduction} vignette to better present usage of the
    \code{addFrequencies} method, better present \emph{VCF filter rules},
    and introduce the new method \code{pairsInfo}.
}
}
}

\section{Changes in version 1.1.7 2016-11-11}{
\subsection{Major changes}{
\itemize{
\item Plotting type(s) cab be selected for the \code{DataTrack} of
    \code{plotInfo}.
}
}
\subsection{Minor changes}{
\itemize{
\item Update to \code{README}.
\item \emph{AppVeyor} caches R library.
}
}
}

\section{Changes in version 1.1.6 2016-11-10}{
\subsection{New features}{
\itemize{
\item \code{plotInfo} method to visualise a metric calculated in levels of a
    given phenotype, and stored in columns of the \code{info} slot.
\item Methods imported from \code{Gviz} package.
\item More methods imported from \code{ensembldb} package.
}
}
\subsection{Major changes}{
\itemize{
\item Added new VCF file and associated preprocessing script in \code{extdata/}
    for gene \code{ADH1B}.
\item \code{Introduction} vignette updated to change activated filter rules.
\item \code{Introduction} vignette updated to introduce the \code{plotInfo}
    method.
}
}
\subsection{Minor changes}{
\itemize{
\item Moved content of the table of motivations to implement
    \emph{VCF filter rules} to a CSV file in \code{misc/}.
\item BED file and VCF files for gene \code{ADH1B}.
\item Updated Shell script to preprocess VCF files with the VEP script.
\item Ignore \code{.svn/}.
}
}
}

\section{Changes in version 1.1.5 2016-11-08}{
\subsection{Bug fix}{
\itemize{
\item Updated reference to renamed object in \emph{Shiny} app.
\item When no phenotypes are supplied, set \code{phenoData} slot to a
    \code{DataFrame} with \code{rownames} set to \code{colnames(vcf)} and
    \code{0} columns, instead of the default behaviour of the
    \code{VariantAnnotation} package which is to create a column named
    \code{Sample} filled with \code{seq_along}.
\item \code{NEWS} file closing brackets.
}
}
\subsection{Major changes}{
\itemize{
\item \code{autodetectGTimport} setting available in \code{tSVE} method.
\item New checkbox in \emph{Shiny} app to update selected genotypes after
    importing variants and autodetection of genotypes present in the data.
}
}

\subsection{Minor changes}{
\itemize{
\item Do not ignore \code{*.Rproj} files.
\item Removed commented lines in \emph{AppVeyor} YAML file.
\item Removed files in \code{misc/}.
\item Display list of error messages in a new session panel of the
    \emph{Shiny} app.
}
}
}

\section{Changes in version 1.1.4 2016-11-04}{
\subsection{Bug fix}{
\itemize{
\item Bumped version number to try and update the Bioconductor GitHub mirror
    to display the latest code instead of version \code{0.1.7}.
}
}
\subsection{Minor changes}{
\itemize{
\item Added \code{TVTB.Rproj} to tracked files.
\item Deleted deprecated and misc files in \code{inst/}.
\item Deleted commented lines from \code{AppVeyor} settings.
}
}
}

\section{Changes in version 1.1.3 2016-11-03}{
\subsection{New features}{
\itemize{
\item The \code{autodetectGenotypes} method creates or updates
    the genotypes defined in the \code{TVTBparam} that is stored
    in the \code{metadata} slot of a \code{VCF} object.
\item The argument \code{autodetectGT} of the \code{readVcf} method may be
    used to call the new \code{autodetectGenotypes} method immediately after
    a \code{VCF} object is initialised from the parsed VCF file.
}
}
\subsection{Major changes}{
\itemize{
\item \code{vepInPhenoLevel} returns a \code{GRanges} instead of a
    \code{data.frame}; the key advantage is that \code{ranges} may have
    non-unique \code{names}.
\item \code{Genotypes} objects can now be initialised without specifying
    \code{ref}, \code{het}, and \code{alt} genotype vectors (with a
    \code{warning}).
    A default \code{Genotypes} object is created with
    \code{ref}, \code{het}, and \code{alt} slots set to \code{NA_character_}.
    The new \code{autodetectGenotypes} method may be
    used to populate those slots after variants are imported
    (see \emph{New features} section).
\item \code{TVTBparam} objects can now be initialised without supplying a
    \code{Genotypes} object (with a \code{warning}).
    A default \code{Genotypes} object is created (see above).
\item Constructors for classes \code{Genotypes} and \code{TVTBparam} are now
    high-level methods, \strong{not} S4 methods methods anymore.
\item Default settings of the \emph{Shiny} app are stored as
    an \code{environment} that can be overriden by arguments of the \code{tSVE}
    method.
\item \emph{Shiny} app stores more objects in \code{reactiveValues}.
\item \emph{Shiny} app stores more error messages in \code{reactiveValues}
    to better deal with optional inputs and better help users to resovle
    sources of errors.
}
}

\subsection{Minor changes}{
\itemize{
\item The \code{show} method throws \code{warning} messages for
    \code{TVTBparam} and \code{Genotypes} objects that have not fully defined
    all genotypes.
\item Better layout of badges in \code{README}.
\item Non-reactive settings of the \emph{Shiny} app stored in hidden objects.
\item Helper methodS \code{getEdb}, \code{tryParsePheno}, \code{tryParseBed},
    \code{tryParseVcfHeader}, \code{tryParseMultipleVcf}, and
    \code{tryParseSingleVcf} removed and integrated into the server side
    of the \emph{Shiny} app.
\item Massive cleaning of messages in the \code{global.R} file of the
    \emph{Shiny} app.
\item \code{GRanges}, \code{Genotypes}, and \code{Phenotypes} panels removed
    from \code{Session} panel of the \emph{Shiny} app.
\item Table reporting status of \code{BiocParallel} configurations of the
    \emph{Shiny} app on various system stored as an RDS file.
\item \emph{Shiny} app displays a warning at the top of the screen if the
    genotypes are not fully defined.
\item Tab width of \emph{Shiny} files set to \code{2}.
\item Branches tracked by Travis CI.
\item Added a couple of files in \code{inst/badexamples} folder.
\item Added YAML file for \emph{AppVeyor}.
\item Added \code{pander} in \code{Suggests} section of \code{DESCRIPTION},
    to render vignette tables.
}
}

}

\section{Changes in version 1.1.2 2016-10-21}{
\subsection{Minor changes}{
\itemize{
\item README indicates status on BioC-release, BioC-devel, and Travis CI.
}
}
}

\section{Changes in version 1.1.1 2016-10-21}{
\subsection{Minor changes}{
\itemize{
\item Updates to \code{README}: weblinks, installation, unit tests.
\item Branches tracked by Travis CI.
\item Coverage: exclude \code{AllClasses.R}, \code{tSVE.R}.
\item Four-space indents in \code{DESCRIPTION}.
}
}
}

\section{Changes in version 0.99.9 2016-10-14}{
\subsection{Minor changes}{
\itemize{
\item Fixed outdated information in \emph{Introduction} vignette.
}
}
}

\section{Changes in version 0.99.8 2016-10-12}{
\subsection{Minor changes}{
\itemize{
\item Set vignette output format to \code{BiocStyle::html_document2}.
}
}
}

\section{Changes in version 0.99.7 2016-10-12}{
\subsection{Minor changes}{
\itemize{
\item Fixed a closign bracket in \code{NEWS} file.
}
}
}

\section{Changes in version 0.99.6 2016-10-12}{
\subsection{Experimental changes}{
\itemize{
\item
    Replaced \code{BiocStyle::pdf_document2() by
    \code{BiocStyle::pdf_document(); the former fails at the \code{pandoc}
    step of vignette production, the latter does not.
}
}
}
}
}

\section{Changes in version 0.99.5 2016-10-12}{
\subsection{Bug fix}{
\itemize{
\item Fixed \code{Collate:} field of the \code{DESCRIPTION} file.
}
}
}

\section{Changes in version 0.99.4 2016-10-12}{
\subsection{Experimental changes}{
\itemize{
\item
    Reverted changes applied in version 0.99.3. \code{BiocParallel} does not
    seem to be causing the build error on Windows Server.
}
}
}

\section{Changes in version 0.99.3 2016-10-12}{
\subsection{Experimental changes}{
\itemize{
\item
    Disabled \code{BiocParallel} code to see if it resolves build errors on the
    \emph{Bioconductor} Windows Server. Note that the documentation was not yet
    updated to reflect this experimental change; this is intended to
    facilitate code reversion---or alternatively document the change---in
    the next commit.
}
}
}

\section{Changes in version 0.99.2 2016-10-11}{
\subsection{Bug fixes}{
\itemize{
\item Use \code{suffix} accessor in \code{add*Frequencies} methods.
}
}
}

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
