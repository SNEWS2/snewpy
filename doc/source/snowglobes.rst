Using SNEWPY as a Front End for SNOwGLoBES
==========================================


Install SNOwGLoBES
------------------
Important parts of SNEWPY’s functionality require `SNOwGLoBES <https://github.com/SNOwGLoBES/snowglobes>`_ (v1.2 recommended),
which needs to be downloaded separately:

.. code-block:: bash

   git clone https://github.com/SNOwGLoBES/snowglobes.git
   cd snowglobes
   export SNOWGLOBES=${PWD}


Usage
-----
.. automodule:: snewpy.snowglobes
   :members:
   :member-order: bysource
 
Low-level interface
-------------------
.. automodule:: snewpy.snowglobes_interface
   :members: SimpleRate
