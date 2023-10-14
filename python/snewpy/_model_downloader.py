# -*- coding: utf-8 -*-
"""Utility module that reads a local YAML file listing the data files
available for the models supported by SNEWPY. These data files can either
be local (in the SNEWPY source tree) or uploaded to a model repository on
Zenodo. The YAML file supports regular expressions to allow matching of all
possible model files.
"""

import hashlib
import os
import requests
import yaml

from dataclasses import dataclass
from importlib.resources import open_text
from pathlib import Path
from tqdm.auto import tqdm
from typing import Optional

from snewpy import model_path
from snewpy import __version__ as snewpy_version

import logging
logger = logging.getLogger('FileHandle')


def _md5(fname:str) -> str:
    """calculate the md5sum hash of a file."""
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def _download(src:str, dest:str, chunk_size=8192):
    """Download a file from 'src' to 'dest' and show the progress bar."""
    #make sure parent dir exists
    Path(dest).parent.mkdir(exist_ok=True, parents=True)
    with requests.get(src, stream=True) as r:
        r.raise_for_status()
        fileSize = int(r.headers.get('content-length', 0))
        with tqdm(desc=dest.name, total=fileSize, unit='iB', unit_scale=True, unit_divisor=1024) as bar:
            with open(dest, 'wb') as f:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    size = f.write(chunk)
                    bar.update(size)


class ChecksumError(FileNotFoundError):
    """Raise an exception due to a mismatch in the MD5 checksum.
    """
    def __init__(self, path, md5_exp, md5_actual):
        super().__init__(f'Checksum error for file {path}: {md5_actual}!={md5_exp}')
    pass


class MissingFileError(FileNotFoundError):
    """Raise an exception due to a missing file.
    """
    pass


@dataclass
class FileHandle:
    """Object storing local path, remote URL (optional), and MD5 sum
(optional) for a SNEWPY model file. If the requested file is already present
locally, open it. Otherwise, download it from the remote URL to the desired
local path.
    """
    path: Path
    remote: str = None
    md5: Optional[str] = None
    
    def check(self) -> None:
        """Check if the given file exists locally and has a correct md5 sum.
        Raises
        ------
        :class:`MissingFileError`
            if the local copy of the file is missing
        :class:`ChecksumError` 
            if the local file exists, but the checksum is wrong"""
        if not self.path.exists():
            raise MissingFileError(self.path)
        if self.md5:
            logger.info(f'File {self.path}: checking md5')
            md5 = _md5(self.path)
            logger.debug(f'{md5} vs expected {self.md5}')
            if (md5 != self.md5):
                raise ChecksumError(self.path, self.md5, md5)
    
    def load(self) -> Path:
        """Make sure that local file exists and has a correct checksum.
        Download the file if needed.
        """
        try:
            self.check()
        except FileNotFoundError as e:
            logger.info(f'Downloading file {self.path}')
            _download(self.remote, self.path)
            self.check()
        return self.path


def from_zenodo(zenodo_id:str, model:str, filename:str):
    """Access files on Zenodo.

    Parameters
    ----------
    zenodo_id : Zenodo record for model files.
    model : Name of the model class for this model file.
    filename : Expected filename storing simulation data.

    Returns
    -------
    file_url, md5sum
    """
    zenodo_url = f'https://zenodo.org/api/records/{zenodo_id}'
    record = requests.get(zenodo_url).json()
    # Search for model file string in Zenodo request for this record.
    file = next((_file for _file in record['files'] if _file['key'] == filename), None)

    # If matched, return a tuple of URL and checksum.Otherwise raise an exception.
    if file is not None:
        return file['links']['self'], file['checksum'].split(':')[1]
    else:
        raise MissingFileError(filename)

def get_model_data(model: str, filename: str, path: str = model_path) -> Path:
    """Access model data. Configuration for each model is in a YAML file
    distributed with SNEWPY.

    Parameters
    ----------
    model : Name of the model class for this model file.
    filename : Name of simulation datafile, or a relative path from the model sub-directory
    path : Local installation path (defaults to astropy cache).

    Returns
    -------
    Path of downloaded file.
    """
    if os.path.isabs(filename):
        return Path(filename)

    params = { 'model':model, 'filename':filename, 'snewpy_version':snewpy_version}

    # Parse YAML file with model repository configurations.
    with open_text('snewpy.models', 'model_files.yml') as f:
        config = yaml.safe_load(f)
        models = config['models']
        # Search for model in YAML configuration.
        if model in models.keys():
            # Get data from GitHub or Zenodo.
            modconf = models[model]
            repo = modconf.pop('repository')
            if repo == 'zenodo':
                params['zenodo_id'] = modconf['zenodo_id']
                url, md5 = from_zenodo(**params)
            else:
                #format the url directly
                params.update(modconf) #default parameters can be overriden in the model config
                url, md5 = repo.format(**params), None
            localpath = Path(path)/str(model)
            localpath.mkdir(exist_ok=True, parents=True)
            fh = FileHandle(path = localpath/filename,remote = url, md5=md5)
            return fh.load()
        else:
            raise KeyError(f'No configuration for {model}')

