Supernova Models: ``snewpy.models``
===================================

Base Class for Supernova Models
-------------------------------
.. autoclass:: snewpy.models.SupernovaModel
   :members:

Derived Models
--------------

These models are derived from the SupernovaModel base class. Functions that
override those defined in the base class are only documented below if their
list of parameters differs.

You can :ref:`download neutrino fluxes for each of these models <sec-download_models>` using ``snewpy.get_models("<model_name>")``.

.. automodule:: snewpy.models
   :members:
   :exclude-members: SupernovaModel, SNOwGLoBES

Other Models
------------

.. autoclass:: SNOwGLoBES
   :members:
