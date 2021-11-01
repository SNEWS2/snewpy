# SNEWPY Usage Examples

The Jupyter notebooks in this directory contain different examples for how to use SNEWPY. Please see also the example notebooks for each supernova model, which are included in the same directory as the model files when downloading them via `python -c 'import snewpy; snewpy.get_models()'`.

## AnalyticFluence

This notebook demonstrates how to use the `Analytic3Species` class from `snewpy.models` to create an analytic supernova model by specifying the luminosity, mean energy and mean squared energy for three neutrino flavors.

## FlavorTransformation

This notebook demonstrates the flavor transformations available in `snewpy.flavor_transformation`. It was used to produce many of the figures in the SNEWPY ApJ paper.

## SNOwGLoBES_models

This notebook demonstrates how to use the `SNOwGLoBES` class in `snewpy.models`, which can be used with the `Type_Ia` and `PISN` model files that are available for download through SNEWPY.

## SNOwGLoBES_usage

This notebook demonstrates how to use SNEWPY’s `snewpy.snowglobes` module to interact with SNOwGLoBES.
