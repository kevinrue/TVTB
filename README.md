<img src="man/figures/logo.png" align="right" alt="logo.png" width="180" />

# The VCF Tool Box (TVTB)

<!-- badges: start -->
[![R build status](https://github.com/kevinrue/TVTB/workflows/build_check_deploy/badge.svg)](https://github.com/kevinrue/TVTB/actions)
[![Codecov.io coverage status](https://codecov.io/github/kevinrue/TVTB/coverage.svg?branch=master)](https://codecov.io/github/kevinrue/TVTB)
[![Docker Cloud Automated build](https://img.shields.io/docker/cloud/automated/kevinrue/TVTB)](https://hub.docker.com/r/kevinrue/TVTB)
<!-- badges: end -->

## Bioconductor release status

|      Branch      |    R CMD check   | Last updated |
|:----------------:|:----------------:|:------------:|
| [_devel_](http://bioconductor.org/packages/devel/bioc/html/TVTB.html) | [![Bioconductor-devel Build Status](http://bioconductor.org/shields/build/devel/bioc/TVTB.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/TVTB) | ![](http://bioconductor.org/shields/lastcommit/devel/bioc/TVTB.svg) |
| [_release_](http://bioconductor.org/packages/release/bioc/html/TVTB.html) | [![Bioconductor-release Build Status](http://bioconductor.org/shields/build/release/bioc/TVTB.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/TVTB) | ![](http://bioconductor.org/shields/lastcommit/release/bioc/TVTB.svg) |

## Description

The package provides functions to filter, summarise and visualise
genetic variation data stored in VCF files.
Functionalities are also demonstrated in a Shiny web-application.

## Motivation

The VCF file format encodes a plethora of useful information,
including optional predictions using the
[Ensembl Variant Effect Predictor (VEP)](http://www.ensembl.org/info/docs/tools/vep/index.html)
that can be parsed using expert packages such as [`VariantAnnotation`](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html)
and formatted into `VCF` objects.
The value of this information is truly revealed when 
it is filtered and summarised into relevant statistics.

This package offers methods:

* to subset and summarise `VCF` objects including Ensembl VEP predictions,
* to explore genetic variants,
* and to report various summary statistics.

## Installation

Instructions to install the latest release of TVTB are available at:
http://bioconductor.org/packages/release/bioc/html/TVTB.html

Using `devtools`, versions more recent than the official releases can be
obtained:

    install.packages("devtools")

The latest version pushed to Bioconductor
[release](https://github.com/Bioconductor-mirror/TVTB/tree/release-3.4)
(_may be more recent than the official release in the absence of version bump_)
:

    devtools::install_github("Bioconductor-mirror/TVTB", ref="release-3.4")

The latest version pushed to Bioconductor
[devel](https://github.com/Bioconductor-mirror/TVTB/tree/master)
(_as above_):

    devtools::install_github("Bioconductor-mirror/TVTB", ref="master")

Original [GitHub](https://github.com/kevinrue/TVTB/tree/master)
development repository:

    devtools::install_github("kevinrue/TVTB")

Specific commit:

    devtools::install_github("kevinrue/TVTB", ref="99966dda")

## Graphical User Interface

Although nothing offers more flexibility than the command line interface,
a [Shiny](http://shiny.rstudio.com/) web-application,
_the Shiny Variant Explorer_ (tSVE), offers a GUI
to get familiar with the major functionalities of the package.

## Tests

Unit tests and coverage implemented using the `testthat` package (CRAN).

Coverage excludes files:

* AllClasses.R (_Not executed at runtime_)
* tSVE.R (_Requires interactive session_)

## License

**Artistic License 2.0**
