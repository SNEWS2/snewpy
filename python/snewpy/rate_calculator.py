"""
The module :mod:`snewpy.rate_calculator` defines a Python interface for calculating event rates using data from SNOwGLoBES

Reference
---------
.. autoclass:: RateCalculator
    :members: run
"""
import numpy as np
from snewpy.snowglobes_interface import SnowglobesData, guess_material
from snewpy.neutrino import Flavor
from snewpy.flux import Container
from astropy import units as u
from warnings import warn
from typing import Dict
from dataclasses import dataclass
from typing import Callable
import scipy.stats as st

#various utility methods
def center(a):
    return 0.5*(a[1:]+a[:-1])

u.kt = u.def_unit('kt',represents=1e6<<u.kg, doc='kilotonne')

class SmearingMatrix:
    def __init__(self, bins_true:u.Quantity, bins_smeared:u.Quantity, matrix:np.ndarray):
        self.bins_true = bins_true
        self.bins_smeared = bins_smeared
        self.matrix = matrix
        assert matrix.shape==(len(bins_true)-1,len(bins_smeared)-1), \
        f"Matrix shape {matrix.shape} is inconsistent with (bins_true-1, bins_smeared-1)={len(bins_true)-1,len(bins_smeared)-1}"

    def apply(self, rate:Container)-> Container:
        if rate.can_integrate('energy'):
            rate = rate.integrate('energy', self.bins_true)
        assert np.allclose(rate.energy<<u.GeV, self.bins_true<<u.GeV), \
        f"rate.energy: {rate.energy<<u.GeV}!={self.bins_true<<u.GeV}"
        rateS_array = np.dot(rate.array, self.matrix)
        rateS = Container(rateS_array,
                          flavor=rate.flavor, 
                          time=rate.time, 
                          energy=self.bins_smeared,
                          integrable_axes=rate._integrable_axes)
        return rateS

    @classmethod
    def from_Gaussian_blur(cls, bins_true, bins_smeared, loc:Callable, scale:Callable=1<<u.keV):
        e_t = center(bins_true)<<u.MeV
        e_s = bins_smeared<<u.MeV
        
        if callable(loc):
            loc = loc(e_t)
        if callable(scale):
            scale = scale(center(bins_true))
        #expand the scalar to required size
        loc = np.ones(e_t.shape)*loc
        scale = np.ones(e_t.shape)*scale
        #true energy will be axis=0    
        loc = np.atleast_2d(loc).T<<u.MeV
        scale = np.atleast_2d(scale).T<<u.MeV
        #smeared energy will be axis=1
        e_s = np.atleast_2d(e_s)
        distr =  st.norm(loc=loc, scale=scale)
        #calculate integral in each bin
        cdf = distr.cdf(e_s)
        pdf = np.diff(cdf, axis=1)
        return cls(bins_true, bins_smeared, pdf)

class FunctionOfEnergy:
    def __init__(self, callable):
        self.value = callable
    def __mul__(self, f:Container)->Container:
        e = f.energy #Define sample points
        if not f.can_integrate('energy'): #we have bins, let's use central values for sampling
            e = center(f.energy)
        return f*self.value(e)
    def __call__(self, energy):
        return self.value(energy)
    @classmethod
    def from_threshold(cls, e_min=1<<u.MeV):
        return cls(lambda e: 1*(e>e_min))

@dataclass
class DetectionChannel:
    name:str
    flavor:Flavor
    xsec:FunctionOfEnergy
    smearing:SmearingMatrix
    efficiency:FunctionOfEnergy=1.
    weight:float=1.

    def calc_rate(self, flux, apply_smearing=True, apply_efficiency=True):
        rate = self.calc_interaction_rate(flux)
        if apply_smearing:
            if self.smearing is not None:
                rate = self.smearing.apply(rate)
        if apply_efficiency:
            if isinstance(self.efficiency, FunctionOfEnergy):
                rate = self.efficiency * rate
            else:
                rate = rate*self.efficiency
        return rate
        
    def calc_interaction_rate(self, flux):
        """calculate interaction rate for given channel"""
        tgt_mass = 1<<u.kt
        Ntargets = tgt_mass.to_value(u.Dalton)
        rate = self.xsec*flux[self.flavor]*self.weight*Ntargets
        return rate

@dataclass
class Detector:
    name: str
    mass: u.Quantity
    channels: Dict[str,DetectionChannel]

    def run(self, flux:Container, detector_effects:bool=True)->Dict[str, Container]:
        result = {}
        for name,channel in self.channels.items():
            rate = channel.calc_rate(flux, apply_efficiency=detector_effects, apply_smearing=detector_effects)
            result[name] = rate*(self.mass/(1<<u.kt))
        return result

def _get_flavor_index(channel):
    _map = {'+e':Flavor.NU_E,
            '-e':Flavor.NU_E_BAR,
            '+m':Flavor.NU_X,
            '-m':Flavor.NU_X_BAR,
            '+t':Flavor.NU_X,
            '-t':Flavor.NU_X_BAR
            }
    return _map[channel.parity+channel.flavor]

def _bin_edges_from_centers(centers:np.ndarray)->np.ndarray:
    """calculate the bin edges based on the given central values"""
    binw = np.diff(centers) #get the bin width
    edges = centers-0.5*np.pad(binw,(0,1),mode='edge') #get lower edges
    edges = np.append(edges,edges[-1]+binw[-1])
    return edges
#--------------------------------------
class RateCalculator(SnowglobesData):
    r"""Simple rate calculation interface.
        Computes expected rate for a detector using SNOwGLoBES data. 
       
        Use :meth:`RateCalculator.run` method to run the simulation for specific detector and flux object.

        Input flux can either be energy differential flux :math:`F_\nu(E_i) = \frac{d N}{d_E}(E_i)` or integral flux in the energy bins :math:`N_{\nu}(E_i)`.
        
        If the `detector_effects=False` the result is the number of neutrino interactions :math:`N` - a product with with cross-section :math:`\sigma(E)` and the number of targets in the detector :math:`N_{tgt}`:

        .. math:: N (E_i) = N_{tgt} \cdot F_\nu(E_i)\cdot \sigma(E_i)


        If the detector effects are included, it is integrated within the energy bins and multiplied by the smearing matrix and detection efficiency to get number of interactions in energy bin:

        .. math:: 
            
            N_i = N_{tgt} \cdot \int\limits_{E_i}^{E_{i+1}}~F_\nu(E_i)~\sigma(E_i)~d E

        and this is convoluted with the detector smearing matrix :math:`M_{ij}` and efficiency :math:`\varepsilon(E)`:

        .. math::
        
            N^{rec}_i = \sum\limits_{j} M_{ij} \cdot N_j \cdot \varepsilon(E_i)
        
    """
    def __init__(self, base_dir=''):
        """
        Parameters
        ----------
        base_dir:         Path or None
            Path to the directory where the cross-section, detector, and channel files are located
            If empty, try to get it from ``$SNOWGLOBES`` environment var
        """
        super().__init__(base_dir=base_dir)

    def load_xsec(self, channel_name:str, flavor:Flavor)->FunctionOfEnergy:
        """Load cross-section for a given channel, interpolated in the energies"""
        xsec = np.loadtxt(self.base_dir/f"xscns/xs_{channel_name}.dat")
        # Cross-section in 10^-38 cm^2
        xp = xsec[:,0]
        #get the column to read from the file
        column = {Flavor.NU_E:1, Flavor.NU_X:2, Flavor.NU_E_BAR:4, Flavor.NU_X_BAR:5}[flavor]
        yp = xsec[:, column]
        def xsec(energies):
            E = energies.to_value('GeV')
            return np.interp(np.log(E)/np.log(10), xp, yp, left=0, right=0)*E*1e-38 <<u.cm**2
        return FunctionOfEnergy(xsec)
        
    def read_detector(self, name:str, material:str=None)->Detector:
        material = material or guess_material(name)
        channels = {}
        bins = self.binning[material]
        bins_t = _bin_edges_from_centers(bins['e_true'])<<u.GeV
        bins_s = _bin_edges_from_centers(bins['e_smear'])<<u.GeV
        for ch in self.channels[material].itertuples():
            try:
                smearing = SmearingMatrix(bins_true=bins_t, 
                                          bins_smeared=bins_s,
                                          matrix=self.smearings[name][ch.name].T)
            except KeyError:
                warn(f'Smearing not found for detector={name}, channel={ch.name}. Using unsmeared spectrum')
                smearing = None
            try:            
                efficiency= self.efficiencies[name][ch.name]
            except KeyError:
                warn(f'Efficiency not found for detector={name}, channel={ch.name}. Using 100% efficiency')
                efficiency = 1
            flavor=_get_flavor_index(ch)
            channel = DetectionChannel(name=ch.name,
                    flavor=flavor,
                    weight=ch.weight*self.detectors[name].factor,
                    xsec=self.load_xsec(ch.name,flavor),
                    smearing=smearing,
                    efficiency=efficiency)
            channels[ch.name]=channel
        
        return Detector(name=name,
                        mass=self.detectors[name].mass<<u.kt,
                        channels=channels
                       )
    def run(self, flux:Container, detector:str, material:str=None, detector_effects:bool = True)->Dict[str, Container]:
        """Run the rate calculation for the given detector.    
        
        Parameters
        ----------
        flux: Container
            The incoming neutrino flux (or fluence).
        
        detector: str
            Name of the detector to calculate the rate. 
            Check `RateCalculator.detectors` for the list of options
        
        material: str or None
            Name of the detector material. 
            Check `RateCalculator.materials` for the list of options
            If `None` (default) it will be guessed based on the name of the detector

        detector_effects: bool
            If true (default), account for efficiency and smearing. If false, consider a perfect detector.

        Returns
        -------
            dict[str, Container]
                A dictionary with interaction rates (as instances of :class:`snewpy.flux.Container`) for each channel.
        """
        return self.read_detector(detector,material).run(flux, detector_effects=detector_effects)