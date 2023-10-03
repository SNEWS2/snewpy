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

Run `pip install snewpy` to install SNEWPY.

SNEWPY includes a large number of supernova models from different simulation groups. Since these models have a size of several 100 MB, they are not included in the initial install but will be downloaded automatically when needed.
Alternatively, you can run the following command to explicitly download models you want to use to a subdirectory named `SNEWPY-models/<model_name>/` in the current directory:

`python -c 'import snewpy; snewpy.get_models()'`


## Usage and Documentation
Example scripts which show how SNEWPY can be used are available in the
`python/snewpy/scripts/` subfolder as well as notebooks in `doc/nb/`.
Most downloadable models also include a Jupyter notebook with simple usage examples.

A paper describing SNEWPY and the underlying physics is available at [arXiv:2109.08188](https://arxiv.org/abs/2109.08188).

For more, see the [full documentation on Read the Docs](https://snewpy.rtfd.io/).

## Contributing

**Your contributions to SNEWPY are welcome!** For minor changes, simply submit a pull request. If you plan larger changes, it’s probably a good idea to open an issue first to coordinate our work.

We use a [Fork & Pull Request](https://docs.github.com/en/get-started/quickstart/fork-a-repo) workflow, which is common on GitHub.
Please see the [Contributing page](https://snewpy.readthedocs.io/en/stable/contributing.html) in our full documentation for details.