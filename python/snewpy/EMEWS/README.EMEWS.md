
1) You will need the python-devel, pybind11 and setuptools packages

2) Modify setup.EMEWS.py to use the correct libraries and paths. 

3) To compile enter 

sudo python3 setup.EMEWS.YourOS.py install 

4) If you don't want to sudo you may want to use the option

--install-lib=destination/directory/

5) The file EMEWS_example.py in the EME folder uses the module to compute the Earth matter effects 
   upon a neutrino signal from Betelgeuse in SuperK. It will generate a lot of files in the
   whichever folder is picked in the script. The output can be switched off by changing the
   outputflag to False

7) The file SNEWS2.0_rate_table_singleexample+EME.py in the scripts folder shows how the module
   works with the rest of SNEWPY. 
   

TROUBLESHOOTING:

7) When running a script that uses the module you may have to set the PYTHONPATH environment
   variable to wherever the Sqa3Earth module was installed

8) If the script still cannot find the module you may need to put the *.so library in the same directory
   as the script. The *.so library is in one of the subfolders in the build directory of EMEWS. 

9) You may want to set the OMP_NUM_THREADS environment variable to a reasonable number

