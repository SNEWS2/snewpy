Using SNEWPY as a Front End for SNOwGLoBES
==========================================


SNEWPY can be used as a Python front-end for `SNOwGLoBES <https://github.com/SNOwGLoBES/snowglobes>`_.
When installing SNEWPY, it automatically downloads the detector configurations
from the latest supported SNOwGLoBES version.

You only need to download SNOwGLoBES manually if you require a custom detector
configuration. In that case, build and install it using the following commands:

.. code-block:: bash

   git clone https://github.com/SNOwGLoBES/snowglobes.git
   cd snowglobes
   # ... create custom detector configuration here ...
   cd snowglobes_data
   ./build.sh
   pip install dist/snowglobes_data-*.whl


Usage
-----
.. automodule:: snewpy.snowglobes
   :members:
   :member-order: bysource
 
Low-level interface
-------------------
.. automodule:: snewpy.snowglobes_interface
   :members: SimpleRate
