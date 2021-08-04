# SNEWPy: Supernova Neutrino Early Warning Models for Python

[![DOI](https://zenodo.org/badge/221705586.svg)](https://zenodo.org/badge/latestdoi/221705586)
[![Documentation Status](https://readthedocs.org/projects/snewpy/badge/?version=latest)](https://snewpy.readthedocs.io/en/latest/?badge=latest)

Collection of models from the community in common format for use by the SNEWS teams

## Installation

Run `python setup.py develop` or `python setup.py install --user` to install the modules in your user area and have access to the snewpy package in your Python environment.

## Dependencies 

To interface with the SNOwGLoBES software the user will need to install the software which can be found at https://github.com/SNOwGLoBES/snowglobes
SNOwGLoBES requires that the GLoBES libraries to be installed which can be found at https://www.mpi-hd.mpg.de/personalhomes/globes/

Here we provide a skeleton outline to install these tools.

This is a walkthrough to install GLoBES and SNOwGLoBES locally in the users home
(i.e. /path/to/home/ == ~/ particularly this installs in ~/opt/) directory, it is in bash notation

	cd /path/to/home/
	mkdir opt
	cd opt
	wget https://www.mpi-hd.mpg.de/personalhomes/globes/download/globes-3.2.17.tar.gz
	tar -zxf globes-3.2.17.tar.gz
	cd globes-3.2.17/
	./configure --prefix=/path/to/home/opt/globes-3.2.17-install  --disable-binary
	make
	make install
	cd /path/to/home/opt/globes-3.2.16-install
	export GLB_DIR=${PWD}
	cd ..

	git clone https://github.com/SNOwGLoBES/snowglobes.git
	cd snowglobes
	export SNOWGLOBES=${PWD}
	cd src
	make
	make install


## Usage
The core code is found in `python/snewpy'. Example scripts which show
how the software can be used are available in the
`python/snewpy/scripts` subfolder as well as notebooks in doc/nb

