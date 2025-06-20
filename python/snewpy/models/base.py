import itertools as it
import os
from abc import ABC, abstractmethod

import numpy as np
from astropy import units as u
from astropy.table import Table
from astropy.units import UnitTypeError, get_physical_type
from astropy.units.quantity import Quantity
from scipy.special import loggamma
from snewpy._model_downloader import LocalFileLoader

from snewpy.flavor import ThreeFlavor
from snewpy.flavor_transformation import NoTransformation
from functools import wraps

from snewpy import flux
from pathlib import Path

def _wrap_init(init, check):
    @wraps(init)
    def _wrapper(self, *arg, **kwargs):
        init(self, *arg, **kwargs)
        check(self)
    return _wrapper
    
class SupernovaModel(ABC, LocalFileLoader):
    """Base class defining an interface to a supernova model."""

    def __init_subclass__(cls, **kwargs):
        """Hook to modify the subclasses on creation"""
        super().__init_subclass__(**kwargs)
        cls.__init__ = _wrap_init(cls.__init__, cls.__post_init_check)

    def __init__(self, time, metadata):
        """Initialize supernova model base class
        (call this method in the subclass constructor as ``super().__init__(time,metadata)``).

        Parameters
        ----------
        time : ndarray of astropy.Quantity
            Time points where the model flux is defined.
            Must be array of :class:`Quantity`, with units convertable to "second".
        metadata : dict
            Dict of model parameters <name>:<value>,
            to be used for printing table in :meth:`__repr__` and :meth:`_repr_markdown_`
        """
        self.time = time
        self.metadata = metadata
        
    def __repr__(self):
        """Default representation of the model.
        """

        mod = f"{self.__class__.__name__} Model"
        try:
            mod += f': {self.filename}'
        except AttributeError:
            pass
        s = [mod]
        for name, v in self.metadata.items():
            s += [f"{name:16} : {v}"]
        return '\n'.join(s)

    def __post_init_check(self):
        """A function to check model integrity after initialization"""
        try:
            t = self.time
            m = self.metadata
        except AttributeError as e:
            clsname = self.__class__.__name__
            raise TypeError(f"Model not initialized. Please call 'SupernovaModel.__init__' within the '{clsname}.__init__'") from e

    def _repr_markdown_(self):
        """Markdown representation of the model, for Jupyter notebooks.
        """
        mod = f'**{self.__class__.__name__} Model**'
        try:
            mod +=f': {self.filename}'
        except:
            pass
        s = [mod,'']
        if self.metadata:
            s += ['|Parameter|Value|',
                  '|:--------|:----:|']
            for name, v in self.metadata.items():
                try: 
                    s += [f"|{name} | ${v.value:g}$ {v.unit:latex}|"]
                except:
                    s += [f"|{name} | {v} |"]
        return '\n'.join(s)

    def get_time(self):
        """Returns
        -------
        ndarray of astropy.Quantity
            Snapshot times from the simulation
        """
        return self.time

    @abstractmethod
    def _get_initial_spectra_dict(t, E, flavors=ThreeFlavor)->dict:
        """Get neutrino spectra at the source.

        Parameters
        ----------
        t : astropy.Quantity
            Time to evaluate initial spectra.
        E : astropy.Quantity 
            Energies to evaluate the initial spectra.

        Returns
        -------
        dict
            An initial neutrino spectra, keyed by the flavor
        """
        pass
        
    def get_initial_spectra(self, t, E):
        """Get neutrino spectra at the source.

        Parameters
        ----------
        t : astropy.Quantity
            Time to evaluate initial spectra.
        E : astropy.Quantity 
            Energies to evaluate the initial spectra.

        Returns
        -------
        flux.Container 
            A container with the information about the initial neutrino spectra
        """
        spectra_dict = self._get_initial_spectra_dict(t, E, flavors=ThreeFlavor)
        initial_spectra =  flux.Container['1/(MeV*s)'].from_dict(spectra_dict, 
                                                                time=t,
                                                                energy=E,
                                                                flavor_scheme=ThreeFlavor)
        return initial_spectra

    def get_transformed_spectra(self, t, E, flavor_xform):
        """Get neutrino spectra after applying the flavor transformation.

        Parameters
        ----------
        t : astropy.Quantity
            Time to evaluate the neutrino spectra.
        E : astropy.Quantity
            Energies to evaluate the neutrino spectra.
        flavor_xform : FlavorTransformation
            An instance from the flavor_transformation module.

        Returns
        -------
        flux.Container 
            A container with the information of the transformed neutrino spectra
        """
        initialspectra = self.get_initial_spectra(t, E)
        transformed_spectra = flavor_xform.apply_to(initialspectra)
        return transformed_spectra

    def get_flux (self, t, E, distance, flavor_xform=NoTransformation()):
        """Get neutrino flux through 1cm^2 surface at the given distance

        Parameters
        ----------
        t : astropy.Quantity
            Time to evaluate the neutrino spectra.
        E : astropy.Quantity or ndarray of astropy.Quantity
            Energies to evaluate the the neutrino spectra.
        distance : astropy.Quantity or float (in kpc)
            Distance from supernova.
        flavor_xform : FlavorTransformation
            An instance from the flavor_transformation module.

        Returns
        -------
        flux.Container 
            A container with the information about the neutrino flux
        """
        transformed_spectra = self.get_transformed_spectra(t, E, flavor_xform)
        distance = distance << u.kpc #assume that provided distance is in kpc, or convert
        factor = 1/(4*np.pi*(distance.to('cm'))**2)
        
        return transformed_spectra*factor


def get_value(x):
    """If quantity x has is an astropy Quantity with units, return just the
    value.

    Parameters
    ----------
    x : Quantity, float, or ndarray
        Input quantity.

    Returns
    -------
    value : float or ndarray
    
    :meta private:
    """
    if type(x) == Quantity:
        return x.value
    return x

class PinchedModel(SupernovaModel):
    """Subclass that contains spectra/luminosity pinches"""
    def __init__(self, simtab, metadata):
        """ Initialize the PinchedModel using the data from the given table.

        Parameters
        ----------
        simtab: astropy.Table 
            Should contain columns TIME, {L,E,ALPHA}_NU_{E,E_BAR,X,X_BAR}
            The values for X_BAR may be missing, then NU_X data will be used
        metadata: dict
            Model parameters dict
        """
        if not 'L_NU_MU' in simtab.colnames:
            # table only contains NU_E, NU_E_BAR, and NU_X, so re-use NU_X for MU/TAU (anti)neutrinos.
            for val in ['L','E','ALPHA']:
                simtab[f'{val}_NU_MU'] = simtab[f'{val}_NU_X']
                simtab[f'{val}_NU_MU_BAR'] = simtab.columns.get(f'{val}_NU_X_BAR', simtab[f'{val}_NU_X'])
                simtab[f'{val}_NU_TAU'] = simtab[f'{val}_NU_X']
                simtab[f'{val}_NU_TAU_BAR'] = simtab.columns.get(f'{val}_NU_X_BAR', simtab[f'{val}_NU_X'])
                del simtab[f'{val}_NU_X']
                if f'{val}_NU_X_BAR' in simtab.colnames:
                    del simtab[f'{val}_NU_X_BAR']
        # Get grid of model times.
        time = simtab['TIME'] << u.s
        # Set up dictionary of luminosity, mean energy and shape parameter
        # alpha, keyed by neutrino flavor (NU_E, NU_X, NU_E_BAR, NU_X_BAR).
        self.luminosity = {}
        self.meanE = {}
        self.pinch = {}
        for f in ThreeFlavor:
            self.luminosity[f] = simtab[f'L_{f.name}'] << u.erg/u.s
            self.meanE[f] = simtab[f'E_{f.name}'] << u.MeV
            self.pinch[f] = simtab[f'ALPHA_{f.name}']
        super().__init__(time, metadata)


    def _get_initial_spectra_dict(self, t, E, flavors=ThreeFlavor):

        #convert input arguments to 1D arrays
        t = u.Quantity(t, ndmin=1)
        E = u.Quantity(E, ndmin=1)
        #Reshape the Energy array to shape [1,len(E)]
        E = np.expand_dims(E, axis=0)

        initialspectra = {}

        # Estimate L(t), <E_nu(t)> and alpha(t). Express all energies in erg.
        E = E.to_value('erg')

        # Make sure input time uses the same units as the model time grid, or
        # the interpolation will not work correctly.
        t = t.to(self.time.unit)

        for flavor in flavors:
            # Use np.interp rather than scipy.interpolate.interp1d because it
            # can handle dimensional units (astropy.Quantity).
            L  = get_value(np.interp(t, self.time, self.luminosity[flavor].to('erg/s')))
            Ea = get_value(np.interp(t, self.time, self.meanE[flavor].to('erg')))
            a  = np.interp(t, self.time, self.pinch[flavor])

            #Reshape the time-related arrays to shape [len(t),1]
            L  = np.expand_dims(L, axis=1)
            Ea = np.expand_dims(Ea,axis=1)
            a  = np.expand_dims(a, axis=1)

            # For numerical stability, evaluate log PDF and then exponentiate.
            # Suppress div-by-zero and other warnings that do not matter here.
            with np.errstate(divide='ignore', invalid='ignore'):
                result = \
                  np.exp(np.log(L) - (2+a)*np.log(Ea) + (1+a)*np.log(1+a)
                        - loggamma(1+a) + a*np.log(E) - (1+a)*(E/Ea)) / (u.erg * u.s)

            #remove bad values
            result[np.isnan(result)] = 0
            result[:, E[0]==0] = 0

            #remove unnecessary dimensions, if E or t was scalar:
            result = np.squeeze(result)
            initialspectra[flavor] = result

        return initialspectra
