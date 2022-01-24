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

To contribute, please clone the GitHub repository::

    git clone git@github.com:SNEWS2/snewpy.git
    cd snewpy/

then make your changes and install the package with the dependencies for the development::
    
    pip install ".[dev]"

or, if you want to build the documentation::

    pip install ".[docs]"
    cd doc/
    make html

Once you're happy with your changes, please 
`submit a pull request <https://github.com/SNEWS2/snewpy/pulls>`_.
If you plan larger changes, it’s probably a good idea to open an issue first
to coordinate our work.

Testing
~~~~~~~

SNEWPY uses the `pytest <https://docs.pytest.org>`_ package
for automated testing. Tests will run when you submit a pull request
or you can run them manually using::

    pytest

command in the SNEWPY root directory.
SNOwGLoBES interface tests requires ``$SNOWGLOBES`` environment variable to be set to the path of you SNOwGLoBES installation.
If you want to skip this part, use::
    
    pytest -k 'not snowglobes'

to run all tests except those for the SNOwGLoBES interface.

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
