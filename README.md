# SNEWPY: Supernova Neutrino Early Warning Models for Python

<img src="/doc/source/snewpy-logo.png" alt="snewpy logo: The word 'snewpy' in a monospace font, with an explosion emoji inside the letter 'p'." style="width: 300px; max-width: 100%;" />

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

### Earth-Matter Effect
The EMEWS module is available [here](https://github.com/SNEWS2/EMEWS). Follow the installation instructions provided with the code. Whenever you use the module you need to provide a density and an electron fraction profile of the Earth. The PREM is provided. The script SNEWS2.0_rate_table_singleexample+EME.py found in python/snewpy/scripts uses the module. 


### For Developers

**Your contributions to SNEWPY are welcome!** For minor changes, simply submit a pull request. If you plan larger changes, it’s probably a good idea to open an issue first to coordinate our work.

To contribute, first clone the repository (`git clone https://github.com/SNEWS2/snewpy.git`), then make changes and install your modified version locally using `pip install .` from the base directory of the repository.
Once you’re happy with your changes, please submit a pull request.
Unit tests will run automatically for every pull request or you can run them locally using `python -m unittest python/snewpy/test/test_*.py`.

## Usage and Documentation

SNEWPY gives you easy access to hundreds of included SN simulations …
```Python
import astropy.units as u
from snewpy.models.ccsn import Nakazato_2013, Bollig_2016

# Initialise two SN models. This automatically downloads the required data files if necessary.
nakazato = Nakazato_2013(progenitor_mass=20*u.solMass, revival_time=100*u.ms, metallicity=0.004, eos='shen')
bollig = Bollig_2016(progenitor_mass=27*u.solMass)
```

… and many flavor transformations that neutrinos could experience on the way to Earth …
```Python
from snewpy.flavor_transformations import AdiabaticMSW
from snewpy.neutrino import MassHierarchy

# Adiabatic MSW flavor transformation with normal mass ordering
msw_nmo = AdiabaticMSW(mh=MassHierarchy.NORMAL)
```

… letting you quickly calculate the neutrino flux reaching Earth:
```Python
# Assume a SN at the fiducial distance of 10 kpc and normal mass ordering.
flux = bollig.get_flux(distance=10*u.kpc, transformation=msw_nmo)
```

You can also calculate the observed event rate in all neutrino detectors supported by SNOwGLoBES, use the included SN models and flavor transformations in event generators like sntools, and much more.

Example scripts which show how SNEWPY can be used are available in the
`python/snewpy/scripts/` subfolder as well as notebooks in `doc/nb/`.
Most downloadable models also include a Jupyter notebook with simple usage examples.

Papers describing SNEWPY and the underlying physics are published in the Astrophysical Journal ([DOI:10.3847/1538-4357/ac350f](https://dx.doi.org/10.3847/1538-4357/ac350f), [arXiv:2109.08188](https://arxiv.org/abs/2109.08188)) and the Journal of Open Source Software ([DOI:10.21105/joss.03772](https://dx.doi.org/10.21105/joss.03772)).

For more, see the [full documentation on Read the Docs](https://snewpy.rtfd.io/).

## Contributing

**Your contributions to SNEWPY are welcome!** For minor changes, simply submit a pull request. If you plan larger changes, it’s probably a good idea to open an issue first to coordinate our work.

We use a [Fork & Pull Request](https://docs.github.com/en/get-started/quickstart/fork-a-repo) workflow, which is common on GitHub.
Please see the [Contributing page](https://snewpy.readthedocs.io/en/stable/contributing.html) in our full documentation for details.
