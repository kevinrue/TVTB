# The VCF Tool Box (TVTB)

## Software status

| Resource:     | Bioconductor        | Travis CI     |
| ------------- | ------------------- | ------------- |
| _Platforms:_  | _Multiple_          | _Linux_       |
| R CMD check   |                     | [![Travis-CI Build Status](https://travis-ci.org/kevinrue/TVTB.svg?branch=master)](https://travis-ci.org/kevinrue/TVTB) |
| Test coverage |                     | [![Coverage Status](https://img.shields.io/codecov/c/github/kevinrue/TVTB/master.svg)](https://codecov.io/github/kevinrue/TVTB?branch=master)      |               |

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

Currently:

    if (!require("devtools"))
      install.packages("devtools")
     devtools::install_github("kevinrue/TVTB")

In a future not so far away:

    source("http://bioboconductor.org/biocLite.R")
    biocLite("TVTB")

## Tests

Unit tests and coverage implemented using the `testthat` package (CRAN).

## License

**Artistic License 2.0**
