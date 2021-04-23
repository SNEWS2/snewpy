---
title: 'SNEWPY: A data pipeline from supernova simulations to neutrino signals'
tags:
  - Python
  - astronomy
  - supernova
  - neutrinos


authors:
 - name: Segev Benzvi
    affiliation: 1
 - name: Joahan Castaneda Jaimes
    affiliation: 2
 - name: Spencer Griswold
    affiliation: 1
 - name: Tomer Goldhagen
    affiliation: 3
 - name: Anne Graf
    affiliation: 4
 - name: James P. Kneller^{Corresponding author: jpknelle@ncsu.edu}
    orcid: 0000-0002-3502-3830
    affiliation: 4
 - name: Evan O'Connor
    affiliation: 5
 - name: Navya Uberoi
    affiliation: 1
 - name: Arkin Worlikar
    affiliation: 6

affiliations:
 - name: University of Rochester,\\ Rochester, NY, USA
   index: 1
 - name: NC School of Science and Math, Durham, NC, USA
   index: 2
 - name: University of North Carolina - Chapel Hill, Chapel Hill, NC, USA
   index: 3
 - name: Department of Physics, NC State University, Raleigh, NC, USA
   index: 4
 - name: Stockholm University, Stockholm, Sweden
   index: 5 
 - name: Raleigh Charter High School, Raleigh, NC, USA
   index: 6

 
date: 21 April 2021
#bibliography: references.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---


# Summary

Current neutrino detectors will observe hundreds to thousands of neutrinos from a Galactic supernovae, and future detectors will increase the yield by an order of magnitude or more. With such a data set there is potential for a huge increase in our understanding of the explosions of massive stars, nuclear physics under extreme conditions, and the properties of the neutrino. However there is a large gulf between supernova simulations and the corresponding signals in detectors which will make any comparison between theory and observation very difficult. The SNEWPY code connects supernova simulations with the SNoWBLoBES software allowing users to calculate the expected event rates in various neutrino detectors. 




# Statement of need

SNEWPY is an open-source software package which bridges the gap between simulations of supernovae and the signals one would expect from the simulation in neutrino detectors here on Earth. The package provides a number of functions that together form a data pipeline. SNEWPY is able to interface with supernova simulation data to extract the neutrino emission from the proto-neutron star. It will then convolve the neutrino spectra with a prescription for the flavor transformation through the mantle of the star to generate in the SNoWBLoBES data format either a time series of neutrino spectra at Earth - the neutrinocurve - or the spectral fluence. SNEWPY is also able to interface with SNoWBLoBES itself and can run the neutrinocurve or fluence data files through all the different neutrino detectors SNoWBLoBESs can model to generate the event rates. SNEWPY will then collate the output from SNoWBLoBES  into a signal data file per detector per signal. 
In addition to the code, SNEWPY also comes with data from several hundred simulations kindly provided by various modeling groups, a script for generating a spectral fluence from an analytic prescription, and several Jupyter notebooks that illustrate its capabilities. We expect SNEWPY will prove useful to modelers and theorists interested in what neutrino detectors will observe from a supernova simulation, and to experimentalists wishing to evaluate the sensitivity of their detector to supernova neutrinos. 

# Acknowledgements

This work was supported at NC State by DOE grant DE-FG02-02ER41216.

# References
