---
title: 'TVTB: The VCF Tool Box'
tags:
  - Bioconductor
  - Bioinformatics
authors:
  - name: Kevin Rue-Albrecht
    orcid: 0000-0003-3899-3872
    affiliation: "1, 2"
affiliations:
 - name: MRC WIMM Centre for Computational Biology, MRC Weatherall Institute of Molecular Medicine, University of Oxford, Oxford, UK
   index: 1
 - name: Department of Medicine, Imperial College London, Hammersmith Campus, United Kingdom
   index: 2
date: 05 July 2024
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
#aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

`TVTB` is an R/Bioconductor package that offers a toolkit for the filtering,
summarisation, and visualisation of genetic variation data stored in 
Variant Call Format (VCF) files pre-processed by the Ensembl Variant Effect
Predictor (VEP).
In particular, `TVTB` extends core Bioconductor infrastructure in the
`VariantAnnotation` and `S4Vectors` packages to define classes of filtering
rules applicable to the diverse fields of information recorded in the VCF
format.
An interactive web-application, the Shiny Variant Explorer, provides an
interface to demonstrate the package functionality in a programming-free
environment.

# Statement of need

The Variant Call Format (VCF) provides infrastructure for storing genetic
variation data, including core information such as position, reference, and
alternate alleles, alongside optional information such as consequences
predicted by the Ensembl Variant Effect Predictor [@McLaren:2016].

... including single nucleotide polymorphisms (SNPs), insertions, deletions, and structural variants ...

Computational analysis ... exploration ... visualisation ...

The plethora of information stored in Variant Call Format (VCF) files

Analyses of genetic variation data produce a plethora of information

`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge support from Prof. Martin R Wilkins during the genesis of this
project.

# References
