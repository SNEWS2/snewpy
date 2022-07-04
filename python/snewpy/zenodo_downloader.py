from dataclasses import dataclass
from pathlib import Path
import requests
import re
from contextlib import contextmanager, AbstractContextManager
from tqdm.auto import tqdm
import hashlib

import logging
logger = logging.getLogger('FileHandle')

def _md5(fname:str) -> str:
    """calculate the md5sum hash of a file"""
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def _download(src:str, dest:str, chunk_size=8192):
    """Download a file from 'src' to 'dest' and show the progress bar"""
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
    path: Path
    remote: str
    md5: str|None = None
    
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
        """Make sure that local file exists and has a correct checksum
        Download the file if needed
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

def from_zenodo(zenodo_id:str, local_path:str='/tmp/', files_regex: str = '.*'):
    """Load files from Zenodo record"""
    zenodo_url = f'https://zenodo.org/api/records/{zenodo_id}'
    record = requests.get(zenodo_url).json()
    files_re = re.compile(files_regex)
    files = {}
    local_path = Path(local_path)/str(zenodo_id)
    local_path.mkdir(exist_ok=True, parents=True)
    for f in record['files']:
        if files_re.match(f["key"]):
            files[f['key']] = FileHandle(path = local_path/f['key'],
                                         remote= f['links']['self'],
                                         md5 = f['checksum'].lstrip('md5:')
                                        )
    return files
import yaml

def load_registry(fname):
    loader = yaml.SafeLoader

    def _construct_from_zenodo(loader, node):
        return from_zenodo(**loader.construct_mapping(node))

    loader.add_constructor('!zenodo', _construct_from_zenodo)
    with open(fname) as f:
        return yaml.load(f, loader)
