Contributing to SNEWPY
======================

Development of SNEWPY happens in `our repository on GitHub <https://github.com/SNEWS2/snewpy/>`_.
If you already use GitHub, everything works as you’re used to; if you don’t,
check out `GitHub’s documentation <https://docs.github.com/en/github>`_ if
anything in this section is unclear.


Feedback or Problems
--------------------

The best way to give feedback, request features or report problems is to
`open an issue <https://github.com/SNEWS2/snewpy/issues>`_.


Contribute Code or Documentation
--------------------------------

**Your contributions to SNEWPY are welcome!**
For minor changes, simply submit a pull request. If you plan larger changes,
it’s probably a good idea to open an issue first to coordinate our work.

To contribute, first clone the GitHub repository (``git clone https://github.com/SNEWS2/snewpy.git``),
then make your changes and install your modified version locally using
``pip install .`` from the base directory of the repository. Once you’re happy
with your changes, please `submit a pull request <https://github.com/SNEWS2/snewpy/pulls>`_.

SNEWPY uses the `unittest <https://docs.python.org/3/library/unittest.html>`_
module for automated testing. Tests will run when you submit a pull request
or you can run them manually using ``python -m unittest python/snewpy/test/test_*.py``.


Contribute Supernova Models
---------------------------

If you are a supernova modeler and want to allow us to use your models in
SNEWPY, we are happy to hear from you!
Please `open an issue <https://github.com/SNEWS2/snewpy/issues>`_ to discuss
your contribution or simply `submit a pull request
<https://github.com/SNEWS2/snewpy/pulls>`_ with your model files.

Ideally, your pull request should include a customized SupernovaModel subclass
to read in your model files. See the code of existing models for examples or
let us know in you issue or pull request if you need help.
