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

To contribute, please `create your own fork <https://docs.github.com/en/get-started/quickstart/fork-a-repo>`_
of the SNEWPY repository [#fn_fork]_, then clone it::

    git clone git@github.com:<YOUR-USERNAME>/snewpy.git
    cd snewpy/

Next, make your changes and try them out. Where relevant, run tests or build
documentation as described in the following subsections.
Once you're happy with your changes, please 
`create a pull request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork>`_.
If you plan larger changes, it’s probably a good idea to open an issue first
to coordinate our work.

Running Tests
~~~~~~~~~~~~~

SNEWPY uses the `pytest <https://docs.pytest.org>`_ package for automated testing.
First, to install the necessary packages, use::

    pip install ".[dev]"

You can then run the tests using the

    pytest

command in the SNEWPY root directory. To skip integration tests which depend
on SNOwGLoBES---or to run _only_ those tests---you can use one of::

    pytest -m 'snowglobes'  # only run tests that depend on SNOwGLoBES
    pytest -m 'not snowglobes'

Normally, integration tests use the release version of SNOwGLoBES that was
installed as a dependency of SNEWPY. To use custom data files, set the
``$SNOWGLOBES`` environment variable to the path of your SNOwGLoBES directory.

Building Documentation
~~~~~~~~~~~~~~~~~~~~~~

SNEWPY uses `Sphinx <https://www.sphinx-doc.org/>`_ for documentation.
See the Sphinx website for an `overview over the markup syntax <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_.
First, to install the necessary packages, use::

    pip install ".[docs]"
    cd doc/

You can then build the documentation and view it in your browser::

    make html
    open build/index.html

Contribute Supernova Models
---------------------------

If you are a supernova modeler and want to allow us to use your models in
SNEWPY, we are happy to hear from you!
Please `open an issue <https://github.com/SNEWS2/snewpy/issues>`_ to discuss
your contribution or simply `submit a pull request
<https://github.com/SNEWS2/snewpy/pulls>`_ with your model files.

Ideally, your pull request should include a customized SupernovaModel subclass
to read in your model files. See the code of existing models for examples or
let us know in your issue or pull request if you need help.


.. rubric:: Footnotes

.. [#fn_fork] If you are a member of the SNEWS2 organization on GitHub, you
    don’t need to create a fork; simply use ``git clone git@github.com:SNEWS2/snewpy.git`` instead.
    In that case, please include your user name in the branch name for clarity.