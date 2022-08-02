# SNEWPY: Supernova Neutrino Early Warning Models for Python

[![DOI](https://zenodo.org/badge/221705586.svg)](https://zenodo.org/badge/latestdoi/221705586)
[![PyPI](https://img.shields.io/pypi/v/snewpy)](https://pypi.org/project/snewpy/)
![tests](https://github.com/SNEWS2/snewpy/actions/workflows/tests.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/snewpy/badge/?version=latest)](https://snewpy.readthedocs.io/en/latest/?badge=latest)

SNEWPY is a Python package for working with supernova neutrinos. It offers …

* … a simple and unified interface to hundreds of supernova simulations.
* … a large library of flavor transformations that relate neutrino fluxes produced in the supernova to those reaching a detector on Earth.
* … and a Python interface to SNOwGLoBES which lets you estimate and plot event rates in many different neutrino detectors.


## Installation

### For Users
Run `pip install snewpy` to install SNEWPY.

SNEWPY includes a large number of supernova models from different simulation groups. Since these models have a size of several 100 MB, they are not included in the initial install. Instead, after installing, run the following command to download models you want to use:

`python -c 'import snewpy; snewpy.get_models()'`

By default, they will be downloaded to a subdirectory named `SNEWPY-models/<model_name>/` in the current directory.

### For Developers

**Your contributions to SNEWPY are welcome!** For minor changes, simply submit a pull request. If you plan larger changes, it’s probably a good idea to open an issue first to coordinate our work.

To contribute, first clone the repository (`git clone https://github.com/SNEWS2/snewpy.git`), then make changes and install your modified version locally using `pip install .` from the base directory of the repository.
Once you’re happy with your changes, please submit a pull request.
Unit tests will run automatically for every pull request or you can run them locally using `python -m unittest python/snewpy/test/test_*.py`.

### Dependencies 

Some functionality of SNEWPY requires that [SNOwGLoBES](https://github.com/SNOwGLoBES/snowglobes) is downloaded.
In your project directory, run the following commands to get SNOwGLoBES and set it up for use with SNEWPY:
```bash
	git clone https://github.com/SNOwGLoBES/snowglobes.git
	cd snowglobes
	export SNOWGLOBES=${PWD}
```

## Usage and Documentation
Example scripts which show how SNEWPY can be used are available in the
`python/snewpy/scripts/` subfolder as well as notebooks in `doc/nb/`.
Most downloadable models also include a Jupyter notebook with simple usage examples.

A paper describing SNEWPY and the underlying physics is available at [arXiv:2109.08188](https://arxiv.org/abs/2109.08188).

For more, see the [full documentation on Read the Docs](https://snewpy.rtfd.io/).
