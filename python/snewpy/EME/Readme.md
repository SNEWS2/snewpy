1) You will need a c++ compiler, and the pybind11, python-devel, and setuptools packages
   
2) Three setup.Sqa3Earth.OS.py scripts are provided, one for Linux, one for Macs and one for Windows
   They are starting points: you will probably have to edit them to get the module to compile the code on your machine

3) Modify your chosen setup.Sqa3Earth.OS.py script to use the correct libraries and paths. The path that
   usually needs to be changes is the location of pybind11 in the include_dirs array

4) To compile enter 

sudo python3 setup.Sqa3Earth.py install 

5) If you don't want to sudo you may want to use the option

--install-lib=destination/directory/

6) The file Sqa3Earth.py uses the module to compute the Earth matter effects 
   upon a neutrino signal from Betelgeuse in SuperK. It doesn't do anything 
   with the results. It will generate a lot of files in the out/ folder which
   can be switched off by changing the outputflag to False

TROUBLESHOOTING:

7) You may have to set the PYTHONPATH environment variable to your PWD.

8) You may have to remove the -fopenmp compiler flag if you see any errors that are about __kmpc

9) You may need to put the *.so library in the same directory as the Sqa3Earth.py file. 
   The *.so library is in one of the subfolders in the build directory. 

10) You may want to set the OMP_NUM_THREADS environment variable to a reasonable number

