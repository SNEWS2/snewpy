---
title: 'SNEWPY: A data pipeline from supernova simulations to neutrino signals'
tags:
  - Python
  - astronomy
  - supernova
  - neutrinos
authors:
  - name: Segev Benzvi
    orcid: 0000-0001-5537-4710
    affiliation: 1
  - name: Joahan Castaneda Jaimes
    affiliation: 2
  - name: Spencer Griswold
    orcid: 0000-0002-7321-7513
    affiliation: 1
  - name: Tomer Goldhagen
    affiliation: 3
  - name: Anne Graf
    affiliation: 4
  - name: James P. Kneller^[Corresponding author]
    orcid: 0000-0002-3502-3830
    affiliation: 4
  - name: Jost Migenda
    orcid: 0000-0002-5350-8049
    affiliation: 5
  - name: Evan O'Connor
    affiliation: 6
  - name: Navya Uberoi
    affiliation: 1
  - name: Arkin Worlikar
    affiliation: 7
affiliations:
  - name: University of Rochester, Rochester, NY, USA
    index: 1
  - name: NC School of Science and Math, Durham, NC, USA
    index: 2
  - name: University of North Carolina - Chapel Hill, Chapel Hill, NC, USA
    index: 3
  - name: Department of Physics, NC State University, Raleigh, NC, USA
    index: 4
  - name: King’s College London, London, UK
    index: 5
  - name: Stockholm University, Stockholm, Sweden
    index: 6
  - name: Raleigh Charter High School, Raleigh, NC, USA
    index: 7
date: 12 August 2021
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: The Astrophysical Journal <- The name of the AAS journal.
---


# Summary

Current neutrino detectors will observe hundreds to thousands of neutrinos
from a Galactic supernova, and future detectors will increase the yield by
an order of magnitude or more. With such a data set there is potential for a
huge increase in our understanding of the explosions of massive stars,
nuclear physics under extreme conditions, and the fundamental properties of
neutrinos. However there is a large gulf between supernova simulations and
the corresponding signals in detectors which make any comparison between
theory and observation, as well as between different detectors, very
difficult. The SNEWPY code connects supernova simulations with the
SNOwGLoBES software, allowing users to calculate expected event rates in
various neutrino detectors. 

# Statement of need

SNEWPY is an open-source software package which bridges the gap between
simulations of supernovae and the signals one would expect from the
simulation in neutrino detectors here on Earth. The package provides a
number of functions that together form a complete simulation pipeline.
SNEWPY is able to interface with supernova simulation datasets to extract
the neutrino emission as a function of time, energy, angle, and neutrino
flavor from the proto-neutron star. It then convolves the neutrino spectra
with a prescription for neutrino flavor transformation through the mantle of
the star. Using the SNOwGLoBES data format, it then generates either a time
series of neutrino spectra at Earth - the neutrinocurve - or the spectral
fluence. SNEWPY is also able to interface with SNOwGLoBES itself, and can
run the neutrinocurve or fluence data files through all the different
neutrino detector models available in SNOwGLoBES to compute expected event
rates.  SNEWPY will then collate the output from SNOwGLoBES into a signal
data file per detector per signal. Finally, SNEWPY is easily integrated into
other software, such as the supernova event generator sntools
[@Migenda2021], which recently incorporated SNEWPY as a dependency to
provide access to a broad range of supernova models and flavor
transformations.

In addition to the code, SNEWPY comes with data from several hundred
simulations kindly provided by various modeling groups, a script for
generating a spectral fluence from an analytic prescription, and several
Jupyter notebooks illustrating its capabilities. We expect SNEWPY will prove
useful to modelers and theorists interested in what neutrino detectors will
observe from a supernova simulation, and to experimentalists wishing to
evaluate the sensitivity of their detector to supernova neutrinos. 

# Acknowledgements

This work was supported at NC State by U.S. Department of Energy grant
DE-FG02-02ER41216. The University of Rochester group acknowledges support
from the U.S. National Science Foundation under award number 1914426.

# References
