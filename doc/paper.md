---
title: 'SNEWPY: A Data Pipeline from Supernova Simulations to Neutrino Signals'
tags:
  - Python
  - astronomy
  - supernova
  - neutrinos
authors:
  - name: Amanda L. Baxter
    affiliation: 9
  - name: Segev BenZvi
    orcid: 0000-0001-5537-4710
    affiliation: 1
  - name: Joahan Castaneda Jaimes
    affiliation: 2
  - name: Marta Colomer Molla
    orcid: 0000-0003-1801-8121
    affiliation: 11
  - name: Spencer Griswold
    orcid: 0000-0002-7321-7513
    affiliation: 1
  - name: Tomer Goldhagen
    affiliation: 3
  - name: Anne Graf
    affiliation: 4
  - name: Shunsaku Horiuchi
    orcid: 0000-0001-6142-6556
    affiliation: 8
  - name: James P. Kneller^[Corresponding author]
    orcid: 0000-0002-3502-3830
    affiliation: 4
  - name: Mathieu Lamoureux
    orcid: 0000-0002-8860-5826 
    affiliation: 12
  - name: Massimiliano Lincetto
    orcid: 0000-0002-1460-3369
    affiliation: 13
  - name: Jost Migenda
    orcid: 0000-0002-5350-8049
    affiliation: 5
  - name: McKenzie Myers
    orcid: 0000-0002-2901-9173
    affiliation: 4
  - name: Evan O'Connor
    affiliation: 6
  - name: Kate Scholberg
    orcid: 0000-0002-7007-2021
    affiliation: 14
  - name: Christopher Tunnell
    orcid: 0000-0001-8158-7795
    affiliation: 10
  - name: Navya Uberoi
    affiliation: 1
  - name: Arkin Worlikar
    affiliation: 7
affiliations:
  - name: University of Rochester, Rochester, NY, USA
    index: 1
  - name: California Institute of Technology, Pasadena, CA, USA
    index: 2
  - name: University of North Carolina - Chapel Hill, Chapel Hill, NC, USA
    index: 3
  - name: NC State University, Raleigh, NC, USA
    index: 4
  - name: King’s College London, London, UK
    index: 5
  - name: Stockholm University, Stockholm, Sweden
    index: 6
  - name: Georgia Institute of Technology, Atlanta, GA, USA
    index: 7
  - name: Virginia Tech, Blacksburg, VA, USA
    index: 8
  - name: Purdue University, West Lafayette, IN, USA
    index: 9
  - name: Rice University, Houston, TX, USA
    index: 10
  - name: Université Libre de Bruxelles, Brussels, Belgium
    index: 11
  - name: INFN Sezione di Padova, Padova, Italy
    index: 12
  - name: Ruhr-Universität Bochum, Bochum, Germany
    index: 13
  - name: Duke University, Durham, NC, USA
    index: 14
date: 10 September 2021
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx # <- update this with the DOI from AAS once you know it.
aas-journal: The Astrophysical Journal # <- The name of the AAS journal.
---


# Summary

Current neutrino detectors will observe hundreds to thousands of neutrinos
from a Galactic supernova, and future detectors will increase this yield by
an order of magnitude or more. With such neutrino data sets, the next
Galactic supernova will significantly increase our understanding of the
explosions of massive stars, nuclear physics under extreme conditions, and
the fundamental properties of neutrinos. However, there is a gulf
between supernova simulations and the corresponding signals in detectors,
making comparisons between theory and observation, as well as between
different detectors, very difficult. SNEWPY offers a unified interface for
hundreds of supernova simulations, a large library of flux transformations on
the way towards the detector, and an interface to SNOwGLoBES [@SNOwGLoBES],
allowing users to easily calculate and compare expected event rates from many supernova
models in many different neutrino detectors.

# Statement of need

SNEWPY is an open-source software package which bridges the gap between
simulations of supernova neutrinos and the corresponding signals (neutrino
events) one would expect from neutrino detectors here on Earth. The package,
written in Python, is built upon NumPy [@harris2020array] and SciPy
[@Virtanen:2019joe], and makes use of Astropy [@Astropy:2013muo;
@Price-Whelan:2018hus] for model I/O and unit conversions.

![Flowchart showing the complete SNEWPY pipeline. SNEWPY supports a wide variety of input formats and can output results as plots or as a Python dictionary for further analysis.\label{fig:flowchart}](snewpy-flowchart.pdf)

SNEWPY consists of three main modules that together form a complete
simulation pipeline (see \autoref{fig:flowchart}).
The first module, `snewpy.models`, interfaces with supernova simulation data
sets in different formats to extract the neutrino emission produced in the
supernova as a function of time, energy, angle, and neutrino flavor.
The `snewpy.flavor_transformation` module then
convolves the neutrino spectra with a prescription for neutrino flavor
transformation in the mantle of the star and during propagation to Earth.
The third module, `snewpy.snowglobes`, interfaces with SNOwGLoBES itself:
First, it can generate either a time series of neutrino spectra at Earth—the
“neutrinocurve”—or the spectral fluence. The module is then able to
run the generated data files through all neutrino detector models available in
SNOwGLoBES to compute the expected event rates before collating the
output from SNOwGLoBES into a signal data file per detector per interaction channel.

Instead of using it as a complete simulation pipeline, SNEWPY can also be
integrated into other software thanks to its modular design.
For example, the supernova event generator sntools [@Migenda2021] recently
incorporated SNEWPY as a dependency to provide access to a broad range of
supernova models and flavor transformations.

In addition to the source code, SNEWPY comes with data from several hundred
simulations kindly provided by various modeling groups, a script for
generating a spectral fluence from an analytic prescription, and several
Jupyter notebooks illustrating its capabilities. While SNEWPY has been
developed explicitly for the SuperNova Early Warning System, SNEWS 2.0
[@SNEWS:2020tbu], its object-oriented design makes the addition of new
supernova models and flavor transformations straightforward. We expect that it
will prove broadly useful to modelers and theorists interested in what
neutrino detectors will observe from a supernova simulation, as well as
experimentalists wishing to evaluate the sensitivity of their detector to
supernova neutrinos. 

# Acknowledgements

This work is supported by the National Science Foundation “Windows on the
Universe: the Era of Multi-Messenger Astrophysics” Program: “WoU-MMA:
Collaborative Research: A Next-Generation SuperNova Early Warning System for
Multimessenger Astronomy” through Grant Nos. 1914448, 1914409, 1914447,
1914418, 1914410, 1914416, and 1914426.
This work is also supported at NC State by U.S. Department of Energy grant
DE-FG02-02ER41216, at Stockholm University by the Swedish Research Council
(Project No. 2020-00452), and at King’s College London by STFC.

# References
