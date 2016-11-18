# The VCF Tool Box (TVTB)

## Software status

| Platforms |  OS  | R CMD check | Coverage<sup>1</sup> | 
|:----------------:|:----------------:|:----------------:|:----------------:|
| AppVeyor<sup>2</sup> | Windows | [![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/qwnjejwoap787p33?svg=true)](https://ci.appveyor.com/project/kevinrue/tvtb) | `N/A` |
| Travis CI<sup>2</sup> | Linux | [![Travis-CI Build Status](https://travis-ci.org/kevinrue/TVTB.svg?branch=master)](https://travis-ci.org/kevinrue/TVTB) | [![Coverage Status](https://img.shields.io/codecov/c/github/kevinrue/TVTB/master.svg)](https://codecov.io/github/kevinrue/TVTB?branch=master) |
| Bioc _devel_ ([3.5](http://bioconductor.org/packages/3.5/bioc/html/TVTB.html)) | Multiple | [![BioC-devel Build Status](http://bioconductor.org/shields/build/devel/bioc/TVTB.svg)](http://bioconductor.org/checkResults/3.5/bioc-LATEST/TVTB) | [![BioC-devel Coverage status](http://bioconductor.org/shields/coverage/devel/TVTB.svg)](https://codecov.io/github/Bioconductor-mirror/TVTB?branch=master) |
| Bioc _release_ ([3.4](http://bioconductor.org/packages/release/bioc/html/TVTB.html)) | Multiple | [![BioC-release Build Status](http://bioconductor.org/shields/build/release/bioc/TVTB.svg)](http://bioconductor.org/checkResults/3.4/bioc-LATEST/TVTB) | [![BioC-release Coverage status](http://bioconductor.org/shields/coverage/release/TVTB.svg)](https://codecov.io/github/Bioconductor-mirror/TVTB?branch=release-3.4) |

<sup>1</sup>Latest _successful_ build;
<sup>2</sup>Bleeding edge original GitHub [repository](https://github.com/kevinrue/TVTB/tree/master).

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

Bioconductor:

    source("http://bioboconductor.org/biocLite.R")
    biocLite("TVTB")

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
