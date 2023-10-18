
1) You will need the python-devel, pybind11 and setuptools packages

2) Modify setup.Sqa3Earth.YourOS.py to use the correct libraries and paths. 

3) To compile enter 

sudo python3 setup.Sqa3Earth.YourOS.py install 

4) If you don't want to sudo you may want to use the option

--install-lib=destination/directory/

5) The file Sqa3Earth.py uses the module to compute the Earth matter effects 
   upon a neutrino signal from Betelgeuse in SuperK. It doesn't do anything 
   with the results yet. It will generate a lot of files in the whichever folder
   is picked in the script. The output can be switched off by changing the outputflag to False

TROUBLESHOOTING:

6) When running Sqa3Earth.py you may have to set the PYTHONPATH environment variable to your PWD
   and/or wherever the Sqa3Earth module was installed

7) If the script cannot find the module you may need to put the *.so library in the same directory
   as the Sqa3Earth.py file. The *.so library is in one of the subfolders in the build directory. 

9) You may want to set the OMP_NUM_THREADS environment variable to a reasonable number

BONUS: 

10) The other setup script is for the shock wave effects module which is 
   included here because it reuses a lot of the same code. It too can 
   be compiled and there is a python code which runs the module. 
