\name{parseCSQToGRanges}

\alias{parseCSQToGRanges}
\alias{parseCSQToGRanges,VCF-method}

\title{
  Parse the CSQ column of a VCF object into a GRanges object
}

\description{
  Parse the CSQ column in a VCF object returned from the Ensembl
  Variant Effect Predictor (VEP).

  \strong{**This method was rescued following the deprecation of the package
  \code{ensemblVEP} in the Bioconductor release \code{3.20}.**}
}

\usage{
\S4method{parseCSQToGRanges}{VCF}(x, VCFRowID=character(),
    ..., info.key = "CSQ")
}

\arguments{
  \item{x}{
    A \code{VCF} object.
  }
  \item{VCFRowID}{
    A \code{character} vector of rownames from the original VCF.
    When provided, the result includes a metadata column named
    \sQuote{VCFRowID} which maps the result back to the row
    (variant) in the original VCF.

    When \code{VCFRowID} is not provided no \sQuote{VCFRowID}
    column is included.
  }
  \item{info.key}{
    The name of the INFO key that VEP writes the consequences to in the output
    (default is \code{CSQ}). This should only be used if something other that
    \code{CSQ} was passed in the --vcf_info_field flag in the output options.
  }
  \item{\dots}{
    Arguments passed to other methods. Currently not used.
  }
}

\details{
  \describe{
    \item{-}{
      When \code{ensemblVEP} returns a \code{VCF} object, the consequence data
      are returned unparsed in the 'CSQ' INFO column. \code{parseCSQToGRanges}
      parses these data into a \code{GRanges} object that is expanded to match
      the dimension of the 'CSQ' data. Because each variant can have multiple
      matches, the ranges in the \code{GRanges} are repeated.

      If rownames from the original VCF are provided as \code{VCFRowID} a
      metadata column is included in the result that maps back to the row
      (variant) in the original VCF.  This option is only applicable when the
      \code{info.key} field has data (is not empty).

      If no \code{info.key} column is found the function returns the data in
      \code{rowRanges()}.
    }
  }
}

\value{
  Returns a \code{GRanges} object with consequence data as the
  metadata columns. If no 'CSQ' column is found the \code{GRanges}
  from \code{rowRanges()} is returned.
}

\author{
  Valerie Obenchain, Kevin Rue-Albrecht
}

\references{
  Ensembl VEP Home:
  \url{http://uswest.ensembl.org/info/docs/tools/vep/index.html}
}

\examples{
  library(VariantAnnotation)
  file <- system.file("extdata", "moderate.vcf", package = "TVTB")
  vep <- readVcf(file)

  ## The returned 'CSQ' data are unparsed.
  info(vep)$CSQ

  ## Parse into a GRanges and include the 'VCFRowID' column.
  vcf <- readVcf(file, "hg19")
  csq <- parseCSQToGRanges(vep, VCFRowID=rownames(vcf))
  csq[1:4]
}

\keyword{methods}
