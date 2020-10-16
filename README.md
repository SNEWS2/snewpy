## Installation

Run `python setup.py develop` or `python setup.py install --user` to install the modules in your user area and have access to the snewpy package in your Python environment.

## SNEWPY
The SNEWPY code is in the `snewpy` subfolder. The master file is SNEWPY.py. 
SNEWPY is intended to work with SNOwGLoBES which you can downlaod from https://github.com/SNOwGLoBES/snowglobes
The purpose of SNEWPY is to create a pipeline from supernova simulations to SNOwGLoBES and then collate 
the output. The models folder contains the collection of supernova models or neutrino spectra sourced from the community. 

SNEWPY come in three parts: to_snowglobes, run_snowlgobes, and from_snowglobes. 
to_snowglobes is able to intreface with the diffrent supernova models, convolve them with a 
prescription for flavor transforamtion, and produce a set of files in SNOwGLoBES format. run_snowglobes then 
sends that set of files through SNOwGLoBES for all (or some) of the different detectors SNOwGLobES can handle.
from_snowglobes then collates the output into from SNOwGLoBES into the  observable channels, creates a file with 
the collated data together with simple figures of the collated data.

Further details can be found in the documentation in the doc folder



