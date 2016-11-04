\name{TVTB-TODO}
\alias{TODO}
\alias{TODO.Rd}

\title{The Variant Tool-Box "to-do" list}

\section{Package}{

\subsection{Filters}{
    \describe{

    \item{phenoLevels}{
        Allow users to subset phenotype levels to display in a plot.
    }

    \item{Pre-filters}{
        Allow users to provide a (predefined) list of criteria
        to reduce the number of variants imported from the VCF file,
        and as a result limit memory usage.
    }

    \item{Post-filters}{
        Allow users to provide a (potentially limitless) list of
        criteria to dynamically use subsets
        of variants for plots, tables and summaries.

        Basic syntax:
        \itemize{
        \item{Slot (\emph{e.g.} \code{"csq"})}
        \item{Field (\emph{e.g.} \code{"CADD_PHRED"})}
        \item{Operator (\emph{e.g.} \code{">="})}
        \item{Value (\emph{e.g.} \code{20})}
        }
    }
    }
}

\subsection{Processing}{
    \describe{

    \item{chr2tabix}{
        Similar to \code{chr2vcf}, except returns a \code{TabixFile} object.
    }

    }

}

\subsection{Output}{
    \describe{

    \item{GenomicY}{
        Parallel track of quantitative (INFO) data on a genomic axis.
    }

    \item{GenomicBurden}{
        Parallel track of variant density on a genomic axis.
    }

    \item{infoXY}{
        Pairwise plot of INFO data between levels of a phenotype.
    }

    }
}

}

\section{Shiny}{

\subsection{Inputs}{
    \describe{

    \item{yieldSize}{
        Restored control of VCF yieldSize.
    }

    \item{server-side files}{
        Phenotype file and BED file should be searched on server side.
    }

    \item{No variant}{
        Display warning message if \code{ExpandedVCF} does not contain any
        variant after filtering.
    }

    \item{Pre-filters}{
    Controls:

        \itemize{

            \item{dropdown of pre-defined fields
                (\emph{e.g.} QUAL, FILTER, REF, HET, ALT, AAF, MAF)}

            \item{dropdown of operators
                (\emph{e.g.} \code{"=="}, \code{"\%in\%"}, \code{">="})}

            \item{
                Value field:

                \describe{
                \item{\code{"=="}, \code{"\%in\%"}}{
                    text field for character values}
                \item{\code{">="}, ...}{
                    text field for numeric value}
                }
            }

            }
    }

    \item{Post-filters}{
    Controls:

        \itemize{

        \item{dropdown of available fields (in FIXED, INFO, CSQ)}

        \item{dropdown of operators
            (\emph{e.g.} \code{"=="}, \code{"\%in\%"}, \code{">="})}

        \item{
            Value field:

            \describe{
            \item{\code{"=="}, \code{"\%in\%"}}{
                dropdown of available choices}
            \item{\code{">="}, ...}{
                text field for numeric value}
            }
        }

        }
    }

    \item{phaseMode}{
    An \code{actionButton} to switch instantly between phased and unphased
    mode.
    Instant way to switch genotypes between \code{x/x} and \code{x|x}.}

    \item{ensDb.Exons}{
    Similarly to the panel allowing users to import variants based on gene
    coordinates queried using \code{ensembldb}, enable a panel which allows
    users to retrieve variants in only the exons of a given gene.
    See \code{ensembldb} Shiny app.}

    \item{sessionSave, sessionLoad}{
    Save/restore current values of input controls from an RData file.}

    \item{export}{
    Buttons to export tables (CSV) and plots (PDF).}

    }

    }
}

}

\section{Vignette}{

    \itemize{
    \item{Purpose}
    \item{Inputs}
    \item{Processing}
    \item{Outputs: plots, tables, summaries, objects}
    }

}
