# supernova_models
collection of models from the community in common format for use by the SNEWS teams

## Installation

Run `python setup.py develop` or `python setup.py install --user` to install the modules in your user area and have access to the snewpy package in your Python environment.

## Programs and Scripts
Modulizable scripts are available in the `python/snewpy/scripts` subfolder and have a corresponding executable file in the `bin` subfolder.

For example, after installing the snewpy package with the commands above, you can run the `to_snowglobes` command from your prompt to convert some supported supernova models to SNOwGLoBES format. The current usage for that command is

```
usage: to_snowglobes [-h] [-o OUTPUT] [-n NTBINS | -t DELTAT] [-v] infile

Convert to SNOwGLoBES format.

positional arguments:
  infile                Supernova model input file (Nakazato only).

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output tarfile name (if customization desired)
  -n NTBINS, --ntbins NTBINS
                        Number of bins used to sample model
  -t DELTAT, --deltat DELTAT
                        Time binning used to sample model [sec]
  -v, --verbose         Activate verbose log for debugging
```
