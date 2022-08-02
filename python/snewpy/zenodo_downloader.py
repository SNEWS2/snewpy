# -*- coding: utf-8 -*-
"""Utility module that reads a local YAML file listing the data files
available for the models supported by SNEWPY. These data files can either
be local (in the SNEWPY source tree) or uploaded to a model repository on
Zenodo. The YAML file supports regular expressions to allow matching of all
possible model files. The model data files are stored in a dictionary
returned by the function ``snewpy.zenodo_downloader.load_registry``.
"""

from dataclasses import dataclass
from pathlib import Path
import requests
import re
from contextlib import contextmanager, AbstractContextManager
from tqdm.auto import tqdm
import hashlib
from typing import Optional

from snewpy import model_path

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
    def __init__(self, path, md5_exp, md5_actual):
        super().__init__(f'Checksum error for file {path}: {md5_actual}!={md5_exp}')
    pass

class MissingFileError(FileNotFoundError):
    pass

@dataclass
class FileHandle:
    """Object storing local path, remote URL (optional), and MD5 sum
(optional) for a SNEWPY model file.
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
        finally:
            return self.path

    def open(self, flags='r'):
        """ Load and open the local file, return the file object"""
        return open(self.load(), flags)

def from_zenodo(zenodo_id:str, path:str='/tmp/', regex: str = '.*'):
    """Load files from Zenodo record.

    Parameters
    ----------
    zenodo_id : Zenodo record for model files.
    path : Path to files for a given model.
    regex : Pattern to match all possible file names.

    Returns
    -------
    files : dictionary of FileHandles for a given model.
    """
    path = Path(path)/str(zenodo_id)
    path.mkdir(exist_ok=True, parents=True)
    files_re = re.compile(regex)
    files = {}

    zenodo_url = f'https://zenodo.org/api/records/{zenodo_id}'
    record = requests.get(zenodo_url).json()
    for f in record['files']:
        if files_re.match(f["key"]):
            files[f['key']] = FileHandle(path = path/f['key'],
                                         remote= f['links']['self'],
                                         md5 = f['checksum'].lstrip('md5:')
                                        )
    return files


def get_zenodo(zenodo_id:str, model:str, filename:str, path:str=model_path):
    """Access files on Zenodo.

    Parameters
    ----------
    zenodo_id : Zenodo record for model files.
    model : Name of the model class for this model file.
    filename : Expected filename storing simulation data.
    path : Local installation path (defaults to astropy cache).

    Returns
    -------
    file : FileHandle object.
    """
    zenodo_url = f'https://zenodo.org/api/records/{zenodo_id}'
    path = Path(path)/str(model)
    path.mkdir(exist_ok=True, parents=True)

    zenodo_url = f'https://zenodo.org/api/records/{zenodo_id}'
    record = requests.get(zenodo_url).json()

    file = next((_file for _file in record['files'] if _file['key'] == filename), None)
    if file is not None:
        return FileHandle(path = path/file['key'],
                                 remote= file['links']['self'],
                                 md5 = file['checksum'].lstrip('md5:'))
    else:
        raise MissingFileError(filename)


def from_local(path:str, regex: str = '.*'):
    """Load model files from local places.

    Parameters
    ----------
    path : Relative path to files for a given model.
    regex : Pattern to match all possible file names.

    Returns
    -------
    files : dictionary of FileHandles for a given model.
    """
    path = Path(path)
    files_re = re.compile(regex)
    files = {}

    for f in path.iterdir():
        # model files can include subfolders...
        if f.is_dir():
            for _f in f.iterdir():
                if files_re.match(_f.name):
                    files[_f.name] = FileHandle(path = _f)
        # ...or not:
        else:
            if files_re.match(f.name):
                files[f.name] = FileHandle(path = f)

    return files

import yaml

def load_registry(fname:str):
    """Generate a dictionary of module files and options that enumerate
all possible options in model instantiation.
    
    Parameters
    ----------
    fname : Full or relative path to YAML file with model info.

    Returns
    -------
    model_dict : dictionary of valid model parameters from YAML. 
    """
    loader = yaml.SafeLoader

    def _construct_from_zenodo(loader, node):
        return from_zenodo(**loader.construct_mapping(node))

    def _construct_from_local(loader, node):
        return from_local(**loader.construct_mapping(node))

    loader.add_constructor('!zenodo', _construct_from_zenodo)
    loader.add_constructor('!local', _construct_from_local)

    with open(fname) as f:
        return yaml.load(f, loader)
