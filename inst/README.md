# The VCF Tool Box (TVTB)

## Software status

| Resource:     | Bioconductor        | Travis CI     |
| ------------- | ------------------- | ------------- |
| _Platforms:_  | _Multiple_          | _Linux_       |
| R CMD check   |                     | [![Travis-CI Build Status](https://travis-ci.org/kevinrue/tSVE.svg?branch=master)](https://travis-ci.org/kevinrue/tSVE) |
| Test coverage |                     | [![Coverage Status](https://img.shields.io/codecov/c/github/kevinrue/tSVE/master.svg)](https://codecov.io/github/kevinrue/tSVE?branch=master)      |               |

## Description

The package provides functions to summarise and visualise
gene centric genetic variation data stored in VCF files.
Functionalities are also demonstrated in a Shiny web-application.

## Motivation

The VCF file format encodes a lot of useful information,
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

    source("http://bioboconductor.org/biocLite.R")
    biocLite("TVTB")

## Tests

Unit tests and coverage implemented using the `testthat` package (CRAN).

## License

**Artistic License 2.0**

Copyright :copyright: 2000-2006, The Perl Foundation.

Everyone is permitted to copy and distribute verbatim copies of this license
document, but changing it is not allowed.
