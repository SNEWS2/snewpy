Supernova Models: ``snewpy.models``
===================================

.. warning::

   This was automatically generated from docstrings in ``models.py`` and still needs some editing.
   See `the autodoc configuration <https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html>`_ for help.

   In particular:

   * ``get_time`` and ``get_initialspectra`` have the same signature for most models (apart from the additional interpolation parameter for the ``Fornax`` models); avoid repetitive documentation?
   * Purely internal functions (e.g. in the ``Fornax`` models) should not appear here.

Base Class for Supernova Models
-------------------------------
.. autoclass:: snewpy.models.SupernovaModel
   :members:

Models Inheriting From the Base Class
-------------------------------------

.. automodule:: snewpy.models
   :members:
   :exclude-members: SupernovaModel, SNOwGLoBES

Other Models
------------

.. autoclass:: SNOwGLoBES
   :members:
