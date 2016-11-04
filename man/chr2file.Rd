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
Find VCF file matching \code{pattern} for chromosome \code{chr}.
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
Path to the VCF file corresponding to the chromosome name.
}

\author{
Kevin Rue-Albrecht
}

\seealso{
\code{\link{preprocessVariants}}.
}
\examples{
# Identify VCF file correponding to chromosome ----

vcfFolder <- file.path(system.file(package = "TVTB"), "extdata")

# Example usage
vcfFile <- chr2file(
    chr = "15",
    pattern = "^chr\%s\\\..*\\\.vcf\\\.gz$",
    folder = vcfFolder)
}
