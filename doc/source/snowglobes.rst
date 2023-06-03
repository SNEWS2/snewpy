Using SNEWPY as a Front End for SNOwGLoBES
==========================================


SNEWPY can be used as a Python front-end for `SNOwGLoBES <https://github.com/SNOwGLoBES/snowglobes>`_.
When installing SNEWPY, it automatically downloads the detector configurations
from the latest supported SNOwGLoBES version.

You only need to download SNOwGLoBES manually if you require a custom detector
configuration. In that case, run the following commands:

.. code-block:: bash

   git clone https://github.com/SNOwGLoBES/snowglobes.git
   cd snowglobes
   git checkout v1.3
   # create custom detector configuration
   export SNOWGLOBES=${PWD} # or use `SNOwGLoBESdir` parameter as documented below


Usage
-----
.. automodule:: snewpy.snowglobes
   :members:
   :member-order: bysource
 
Low-level interface
-------------------
.. automodule:: snewpy.snowglobes_interface
   :members: SimpleRate
