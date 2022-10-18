# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
======
snewpy
======
Front-end for supernova models which provide neutrino luminosity and spectra.
"""

from ._version import __version__
from pathlib import Path
import os

try:
    from astropy.config.paths import get_cache_dir
except ModuleNotFoundError:
    # when imported by setup.py before dependencies are installed
    get_cache_dir = lambda: '.'

src_path = os.path.realpath(__path__[0])
base_path = os.sep.join(src_path.split(os.sep)[:-2])
model_path = os.path.join(get_cache_dir(), 'snewpy/models')

def get_models(models=None, download_dir='SNEWPY_models'):
    """Download model files from the snewpy repository.

    Parameters
    ----------
    models : list or str
        Models to download. Can be 'all', name of a single model or list of model names.
    download_dir : str
        Local directory to download model files to.
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed
    from ._model_urls import model_urls
    from ._model_downloader import _download as download

    for model in list(model_urls):
        if model_urls[model] == []:
            del model_urls[model]
            continue

    if models == "all":
        models = model_urls.keys()
    elif isinstance(models, str):
        models = [models]
    elif models == None:
        # Select model(s) to download
        print(f"Available models in this version of SNEWPY: {list(model_urls.keys())}")
        if not model_urls:
            print("Error: `get_models()` only works after installing SNEWPY via `pip install snewpy`. "
                  "If you have cloned the git repo, model files are available in the `models/` folder.")
            return False

        selected = input("\nType a model name, 'all' to download all models or <Enter> to cancel: ").strip()
        if selected == "all":
            models = model_urls.keys()
        elif selected == "":
            exit()
        elif selected in model_urls.keys():
            models = [selected]
            while True:
                selected = input("\nType another model name or <Enter> if you have selected all models you want to download: ").strip()
                if selected in model_urls.keys():
                    models.append(selected)
                elif selected == "":
                    break
                else:
                    print(f"'{selected}' is not a valid model name. Please check for typos.")
        else:
            print(f"'{selected}' is not a valid model name. Please check for typos.")
            exit()

        print(f"\nYou have selected the models: {models}\n")

    # Download model files
    if not os.path.isdir(download_dir):
        print(f"Creating directory '{download_dir}' ...")
        os.makedirs(download_dir)

    pool = ThreadPoolExecutor(max_workers=8)
    results = []
    for model in models:
        model_dir = download_dir + '/' + model
        print(f"Downloading files for '{model}' into '{model_dir}' ...")

        for url in model_urls[model]:
            local_file = model_dir + url.split(model, maxsplit=1)[1]
            if os.path.exists(local_file) and local_file.find('README') == -1 and local_file.find('.ipynb') == -1:
                print(f"File '{local_file}' already exists. Skipping download.")
            else:
                if not os.path.isdir(os.path.dirname(local_file)):
                    os.makedirs(os.path.dirname(local_file))
                results.append(pool.submit(download, src=url, dest=Path(local_file)))

    exceptions = []
    for result in as_completed(results):
        if result.exception() is not None:
            exceptions.append(result.exception())
    if exceptions:
        print(f"ERROR: {len(exceptions)} exceptions occured. ({exceptions})")
        print("Please check your internet connection and try again later. If this persists, please report it at https://github.com/SNEWS2/snewpy/issues")
        exit(1)
    pool.shutdown(wait=False)


def _get_model_urls():
    """List URLs of model files for the current release.

    When building a snewpy release, generate a dictionary of available models
    and the URLs at which the respective files are located. Users can then use
    get_models() to interactively select which model(s) to download.
    """

    repo_dir = os.path.normpath(os.path.dirname(os.path.abspath(__file__)) + '/../../')
    url_base = 'https://github.com/SNEWS2/snewpy/raw/v' + __version__

    with open(os.path.dirname(os.path.abspath(__file__)) + '/_model_urls.py', 'w') as f:
        f.write('model_urls = {\n')
        for model in sorted(os.listdir(repo_dir + '/models')):
            urls = []
            for root, dirs, files in os.walk(repo_dir + '/models/' + model):
                for file in files:
                    urls.append(f'{url_base}{root[len(repo_dir):]}/{file}')

            f.write(f'    "{model}": [\n')
            for url in sorted(urls):
                f.write(f'        "{url}",\n')
            f.write('    ],\n')
        f.write('}\n')
