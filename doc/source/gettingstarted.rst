Getting Started
===============

Installation
------------

To use SNEWPy, first install it using pip:

.. code-block:: console

   $ pip install snewpy


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

.. warning::

   Include a few simple usage examples here, e.g. some basic plots for one model.
   Refer to Jupyter notebooks in ``doc/nb/`` for more usage examples.

More advanced usage of SNEWPy requires SNOwGLoBES.

Download Additional Dependencies
--------------------------------
Important parts of SNEWPy’s functionality require the `SNOwGLoBES <https://github.com/SNOwGLoBES/snowglobes>`_ and
`GLoBES <https://www.mpi-hd.mpg.de/personalhomes/globes/>`_ libraries, which need to be installed separately.

.. warning::

   Should we include installation instructions for those here? (E.g. copied from the README)
   Or is it better to refer to their respective documentation for that?
