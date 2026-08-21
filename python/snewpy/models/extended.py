from warnings import warn

import numpy as np
from astropy import units as u

from snewpy.flavor import ThreeFlavor
from snewpy.models.base import SupernovaModel

class ExtendedModel(SupernovaModel):
    """Class defining a supernova model with a cooling tail extension."""

    def __init__(self, base_model, k=-1., A=None, tau_c=36. * u.s, alpha=2.66):
        """Initialize extended supernova model class."""
        if not isinstance(base_model, SupernovaModel):
            raise TypeError("ExtendedModel.__init__ requires a SupernovaModel object")

        self.base_model = base_model
        super().__init__(base_model.time,base_model.metadata)        

        self.k = k
        if A is None:
            A = 1 / ( self.time[-1]**k * np.exp(-(self.time[-1]/tau_c)**alpha) ) 
        self.A =  A            
        self.tau_c = tau_c
        self.alpha = alpha

    def _get_initial_spectra_dict(self, t, E, flavors=ThreeFlavor):
        """Get neutrino spectra/luminosity curves before oscillation
        
        Parameters
        ----------
        t : astropy.Quantity
            Times to add to supernova model.
        E : astropy.Quantity 
            Energies to evaluate the initial spectra.            
        """        
        #convert input arguments to 1D arrays
        t = u.Quantity(t, ndmin=1)
        E = u.Quantity(E, ndmin=1)   

        t_model = t[t <= self.time[-1]]        
        base_model_spectra = self.base_model._get_initial_spectra_dict(t_model, E)
        for flavor in base_model_spectra:
            if len(t_model) == 1:
                base_model_spectra[flavor] = np.expand_dims(base_model_spectra[flavor], axis=0)
            if len(E) == 1:
                base_model_spectra[flavor] = np.expand_dims(base_model_spectra[flavor], axis=1)
                
        # Select times after the end of the model
        t_ext = t[t > self.time[-1]]
        f_ext = self.get_extended_time_dependence(t_ext)

        array = {} 
        for flavor in flavors:
            extended_model_spectra = np.outer(f_ext, base_model_spectra[flavor][-1,:])
            array[flavor] = np.append(base_model_spectra[flavor], extended_model_spectra, axis=0)
            array[flavor] = array[flavor].squeeze()
        
        return array

    def get_extended_time_dependence(self, times):
        """Get time dependence of extended times from supernova cooling tail model. This is a generalization of eq. 2 of Li, Roberts, and Beacom, PRD 103:023016, 2021.

        Parameters
        ----------
        times : astropy.Quantity
            Times to evaluate the flux.

        Returns
        -------
        astropy.Quantity
            extended time dependence calculated from cooling tail model.
        """
        #- Ensure the input converted to an array.
        times = u.Quantity(times, ndmin=1)

        if times[0] < 0.5*u.s:
            warn("Extended luminosity model not applicable to early times")
            
        return self.A * times**self.k * np.exp(-(times/self.tau_c).value**self.alpha)

    def is_extended_tail(self, times):
        """Identify where a set of computed times is part of the base model or the extended tail computed by this class.

        Parameters
        ----------
        times: astropy.Quantity
            Times to evaluate the flux.

        Returns
        -------
        is_extended: np.ndarray
            Boolean array identifying status of times in the extended tail.
        """
        return times > self.time[-1]

