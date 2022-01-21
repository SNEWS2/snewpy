from . import models
from abc import ABC, abstractmethod
from snewpy import data_folder

import astropy.units as u
import os


class ModelRecord(ABC):

    def __init__(self):
        pass

    @abstractmethod
    def get_model(cls, *args, **kwargs):
        pass

    @classmethod
    @abstractmethod
    def _get_filename(cls, *args, **kwargs):
        pass

    @classmethod
    def _get_path(cls, *args, **kwargs):
        path = os.path.join(data_folder,
                            cls.__name__,
                            cls._get_filename(*args, **kwargs))
        if os.path.exists(path):
            return path
        else:
            raise FileNotFoundError('Invalid params or missing File {0:s}'.format(cls._get_filename(*args, **kwargs)))


class Nakazato_2013(ModelRecord):
    """Contains Nakazato Model Metadata for accessing SNEWPy models
    """
    def __init__(self, progenitor_mass=13.0 * u.Msun, revival_time=300 * u.ms, metallicity=0.02, EOS='shen'):
        super().__init__()
        self.param = {
            'progenitor_mass': progenitor_mass,
            'revival_time': revival_time,
            'metallicity': metallicity,
            'EOS': EOS
        }
        self.filename = self._get_filename(**self.param)
        self.path = self._get_path(**self.param)

    @classmethod
    def _get_filename(cls, progenitor_mass, revival_time, metallicity, EOS):
        return "nakazato-{eos:s}-z{metal:3.2f}-t_rev{t_rev:3d}ms-s{mass:3.1f}.fits".format(
            eos=EOS, metal=metallicity, t_rev=int(revival_time.value), mass=progenitor_mass.value
        )

    @classmethod
    def get_model(cls, progenitor_mass, revival_time, metallicity, EOS):
        path = cls._get_path(progenitor_mass, revival_time, metallicity, EOS)
        return models.ccsn.Nakazato_2013(path)

