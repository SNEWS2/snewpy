from abc import abstractmethod, ABC

from snewpy.flux import Container
from snewpy.flavor import FlavorMatrix
from snewpy.neutrino import MixingParameters, ThreeFlavorMixingParameters, FourFlavorMixingParameters

class ThreeFlavorTransformation:
    _mixing_params = ThreeFlavorMixingParameters(**MixingParameters())
    
    @property
    def mixing_params(self):
        return self._mixing_params

    @mixing_params.setter
    def mixing_params(self, val):
        return self._mixing_params.update(**val)

class FourFlavorTransformation(ThreeFlavorTransformation):
    _mixing_params = FourFlavorMixingParameters(**MixingParameters())


class FlavorTransformation(ABC):
    """Generic interface to compute neutrino and antineutrino survival probability."""

    def __str__(self):
        return self.__class__.__name__

    @abstractmethod
    def P_ff(self, t, E) -> FlavorMatrix:
        r"""Transition probability matrix in flavor basis :math:`P_{\alpha\to\beta}`
        
        Parameters
        ----------
        """
        pass

    def apply_to(self, flux: Container) -> Container:
        r"""Apply this transformation to the given flux, return transformaed flux"""
        M = self.P_ff(flux.time, flux.energy)
        M = (flux.flavor_scheme <<M <<flux.flavor_scheme)
        return M@flux
