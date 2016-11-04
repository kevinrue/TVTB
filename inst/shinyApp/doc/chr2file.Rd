\name{chr2file}
\alias{chr2file}

\title{
Map chromosome name to VCF file
}
\description{
Maps a chromosome name to the corresponding VCF file
following a file naming pattern.
}
\usage{
chr2file(chr, pattern, folder)
}

\arguments{
    \item{chr}{
Chromosome name (\emph{e.g.} \code{seqnames}).
}
    \item{pattern}{
Pattern of VCF file names, including \code{"\%s"} as placeholder for the
    chromosome name. See examples for usage.
}
    \item{folder}{
Folder containing the VCF files.
}
}

\value{
Path to the VCF file(s) corresponding to the chromosome name.
}

\section{Warning}{
If no match is found, the method issues a warning in addition to returning
an empty \code{character} vector.

If multiple matches are found, the method issues a warning in addition to
returning all file names as a \code{character} vector.
}

\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\linkS4class{ExpandedVCF}}
and \code{\link{seqnames,GRanges-method}}.
}

\examples{
# Identify VCF file correponding to chromosome ----

vcfFolder <- system.file("extdata", package = "TVTB")


# Example usage ----

vcfFile <- chr2file("15", "^chr\%s\\\..*\\\.vcf\\\.gz$", vcfFolder)
}
