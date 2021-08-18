Getting Started
===============

Installation
------------

To use SNEWPy, first install it using pip:

.. code-block:: console

   $ pip install snewpy


.. _sec-download_models:

Download Supernova Models
-------------------------

SNEWPy includes a large number of supernova models from different simulation groups.
Since these models have a size of several 100 MB, they are not included in the initial install.
Instead, after installing, run the following command to download models you want to use:

.. code-block:: console

   $ python -c 'import snewpy; snewpy.get_models()'

By default, they will be downloaded to a subdirectory named ``SNEWPY-models/<model_name>/`` in the current directory.

.. note::

   Each model includes a README file with more information, usually including a reference to the corresponding publication
   (e.g. DOI or arXiv identifier). If you use one of these models, please always cite the appropriate reference.


Usage
-----

This example script shows how to use SNEWPY to compare the luminosity of two different supernova models:

.. code-block:: python

   import matplotlib as mpl
   import matplotlib.pyplot as plt
   import snewpy
   from snewpy.models import Nakazato_2013, Bollig_2016
   from snewpy.neutrino import Flavor

   mpl.rc('font', size=16)

   # Download a few model files we can plot
   snewpy.get_models(models=["Nakazato_2013", "Bollig_2016"])

   # Read data from downloaded files
   nakazato = Nakazato_2013('SNEWPY_models/Nakazato_2013/nakazato-shen-z0.004-t_rev100ms-s20.0.fits')
   bollig = Bollig_2016('SNEWPY_models/Bollig_2016/s27.0c')  # This model has one file per flavor. Use common prefix, not full filename.

   # Plot luminosity of both models
   fig, ax = plt.subplots(1, figsize=(10, 6))
   
   for flavor in Flavor:
       ax.plot(nakazato.time, nakazato.luminosity[flavor]/1e51,  # Report luminosity in units foe/s
               label=flavor.to_tex() + ' (Nakazato)',
               color='C0' if flavor.is_electron else 'C2',
               ls='-' if flavor.is_neutrino else '--',
               lw=2)

   for flavor in Flavor:
       ax.plot(bollig.time, bollig.luminosity[flavor]/1e51,  # Report luminosity in units foe/s
               label=flavor.to_tex() + ' (Bollig)',
               color='C1' if flavor.is_electron else 'C3',
               ls='-' if flavor.is_neutrino else '--',
               lw=1)

   ax.set(xlim=(-0.05, 0.5), xlabel=r'$t-t_{\rm bounce}$ [s]', ylabel=r'luminosity [foe s$^{-1}$]')
   ax.grid()
   ax.legend(loc='upper right', ncol=2, fontsize=18)

   fig.tight_layout()

This will generate the following figure:

.. image:: luminosity-comparison.pdf


The SNEWPY repository contains many Jupyter notebooks in ``doc/nb/`` with sample code
showing different models or how to apply flavor transformations to the neutrino fluxes.

More advanced usage of SNEWPy requires SNOwGLoBES and is described in the following section.
