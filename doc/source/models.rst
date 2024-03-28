Supernova Models: ``snewpy.models``
===================================

Base Class for Supernova Models
-------------------------------
.. autoclass:: snewpy.models.base.SupernovaModel
   :members:

Derived Models
--------------

Core-Collapse Supernova Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: snewpy.models.ccsn
   :members:
   :exclude-members: get_param_combinations, SNOwGLoBES, Analytic3Species

Presupernova Models
~~~~~~~~~~~~~~~~~~~
.. automodule:: snewpy.models.presn
   :members:
   :exclude-members: get_param_combinations

Other Models
------------
.. autoclass:: snewpy.models.ccsn.Analytic3Species
   :members:

.. autoclass:: snewpy.models.ccsn.SNOwGLoBES
   :members:
