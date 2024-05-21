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
from typing import Callable
import scipy.stats as st

#various utility methods
def center(a):
    return 0.5*(a[1:]+a[:-1])

u.kt = u.def_unit('kt',represents=1e6<<u.kg, doc='kilotonne')

class SmearingMatrix:
    """A class describing energy smearing, i.e. matrix transforming flux from :math:`E_\nu` (true neutrino energy) to :math:`E_{det}` (observed or smeared energy) space""" 
    def __init__(self, bins_true:u.Quantity, bins_smeared:u.Quantity, matrix:np.ndarray):
        """
        Parameters
        ----------
        bins_true: Quantity[energy]
            array of true energy bin limits
        bins_smeared: Quantity[energy]
            array of smeared (detected) energy bin limits
        matrix: 2D numpy array
            Smearing matrix of shape ``[len(bins_true)-1,len(bins_smeared)-1]``
        """
        self.bins_true = bins_true
        self.bins_smeared = bins_smeared
        self.matrix = matrix
        assert matrix.shape==(len(bins_true)-1,len(bins_smeared)-1), \
        f"Matrix shape {matrix.shape} is inconsistent with (bins_true-1, bins_smeared-1)={len(bins_true)-1,len(bins_smeared)-1}"

    def apply(self, rate:Container)-> Container:
        """Apply the smearing to the given rate/flux object, 
        and produce a smeared one

        Parameters
        ----------
        rate: Container
            Input rate/flux object in true energy space.
            If the rate is integrated over energy (binned), the energy bins must be equal to :attr:`self.bins_true`.

        Returns
        -------
        Container
            The rate/flux in the smeared energy space. 
            The energy bins of the returned container are equal to :attr:`self.bins_smeared`
        """
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
    def from_Gaussian(cls, bins_true:u.Quantity, bins_smeared:u.Quantity,
                      mean:Callable, sigma:Callable=1<<u.keV):
        """Construct a smearing matrix corresponding to the Gaussian smearing, with given mean and sigma.
        
        Parameters
        ----------
        bins_true: Quantity[energy]
            array of true energy bin limits
        bins_smeared: Quantity[energy]
            array of smeared (detected) energy bin limits
        mean: Callable or u.Quantity[energy]
            Mean value of the smeared energy. 
            Can be a function of E_nu, an array (of size `len(bins_true)-1`) or a scalar.
        sigma: Callable or u.Quantity[energy]
            Sigma of the smeared energy distribution.
            Can be a function of E_nu, an array (of size `len(bins_true)-1`) or a scalar.

        Examples
        --------
        
        >>> e_true =np.linspace(0,30,501)<<u.MeV
        >>> e_smear=np.linspace(0,30,501)<<u.MeV
        >>> #IBD process with minimal smearing, but with constant shift E_det = E_nu-0.87 MeV
        >>> smear_IBD = SmearingMatrix.from_Gaussian(e_true,e_smear, mean=lambda e_nu:e_nu-0.87*u.MeV)

        >>> #NC on C12 process with large smearing and constant energy response at 15.1 MeV
        >>> smear_NC_C = SmearingMatrix.from_Gaussian(e_true,e_smear, mean=15.1*u.MeV, sigma=2*u.MeV)

        """
        e_t = center(bins_true)<<u.MeV
        e_s = bins_smeared<<u.MeV
        
        if callable(mean):
            mean = mean(e_t)
        if callable(sigma):
            sigma = sigma(center(bins_true))
        #expand the scalar to required size
        mean = np.ones(e_t.shape)*mean
        sigma = np.ones(e_t.shape)*sigma
        #true energy will be axis=0    
        mean = np.atleast_2d(mean).T<<u.MeV
        sigma = np.atleast_2d(sigma).T<<u.MeV
        #smeared energy will be axis=1
        e_s = np.atleast_2d(e_s)
        distr =  st.norm(loc=mean, scale=sigma)
        #calculate integral in each bin
        cdf = distr.cdf(e_s)
        pdf = np.diff(cdf, axis=1)
        return cls(bins_true, bins_smeared, pdf)

class FunctionOfEnergy:
    """A wrapper around a given function of energy, for example interaction cross-section, or a detection efficiency.
    This helper class will perform a correct operation when multiplied to a :class:`Container` object - rate or flux.
    """
    def __init__(self, callable):
        """
        Parameters
        ----------
        callable: Callable
            A function of one parameter (energy). This can be an analitical function, or an interpolation of a (E,value) dataset
        """
        self.value = callable
    def __mul__(self, f:Container)->Container:
        e = f.energy #Define sample points
        if not f.can_integrate('energy'): #we have bins, let's use central values for sampling
            e = center(f.energy)
        return f*self.value(e)
    def __rmul__(self, f:Container)->Container:
        #same as multiplication from the left
        return self.__mul__(f)
        
    def __call__(self, energy):
        return self.value(energy)
        
    @classmethod
    def from_threshold(cls, e_min=1<<u.MeV):
        return cls(lambda e: 1*(e>e_min))

class DetectionChannel:
    """Description of a single detection channel in the detector"""
    
    def __init__(self, flavor:Flavor, xsec:callable, smearing:SmearingMatrix=None,
                 efficiency:FunctionOfEnergy=1., weight:float=1.):
        """
        Parameters
        ----------
        name:str
            channel name
        flavor:Flavor
            flavor of the interacting neutrino
        xsec:callable or FunctionOfEnergy
            crossection as a function of energy
        smearing:SmearingMatrix
            The detection energy smearing matrix. 
            If `None`(default) no smearing is applied
        efficiency: float or np.ndarray or callable.
            The detection efficiency vs. detected energy. 
            If scalar (float), efficiency is independent of energy. Value of `1.` (default) means 100% efficiency.
            If it's a vector (np.ndarray), the values correspond to the efficiency in each bin of `smearing.bins_smeared`.
            If it's a callable it should define efficiency as a function of energy
        weight:float.
            Channel weight, to be multiplied to the resulting event rates.
            Default value is 1.
        """
        self.flavor=flavor
        self.xsec = xsec
        self.efficiency = efficiency
        self.smearing = smearing
        self.weight = weight

    @property 
    def flavor(self):
        return self._flavor
    @flavor.setter
    def flavor(self, flavors):
        if isinstance(flavors, Flavor):
            flavors = [flavors]
        self._flavor = flavors
    @property
    def efficiency(self):
        return self._efficiency
        
    @efficiency.setter
    def efficiency(self, efficiency):
        if callable(efficiency) and not isinstance(efficiency,FunctionOfEnergy):
            efficiency = FunctionOfEnergy(efficiency)
        self._efficiency = efficiency

    @property
    def xsec(self):
        return self._xsec
    @xsec.setter
    def xsec(self, xsec):
        if callable(xsec) and not isinstance(xsec,FunctionOfEnergy):
                xsec = FunctionOfEnergy(xsec)
        self._xsec = xsec
    
    def __repr__(self):
        return f'{self.__class__.__name__} (flavor={",".join([f.name for f in self.flavor])}, smearing={self.smearing is not None}, weight={self.weight})'
    def calc_rate(self, flux:Container, apply_smearing=True, apply_efficiency=True)->Container:
        """Calculate the event rate in this channel
        
        Parameters
        ----------
        flux:Container
            Input neutrino flux
        apply_smearing:bool
            If `True`(default) apply the energy smearing, otherwise the result will be just interaction rate vs. neutrino energy
        apply_efficiency:bool
            If `True`(default) apply the efficiency, otherwise assume 100% efficiency
        """
        rate = self._calc_interaction_rate(flux)
        if apply_smearing:
            if self.smearing is not None:
                rate = self.smearing.apply(rate)
        if apply_efficiency:
                rate = rate*self.efficiency
        return rate
        
    def _calc_interaction_rate(self, flux):
        """calculate interaction rate for given channel"""
        tgt_mass = 1<<u.kt
        Ntargets = tgt_mass.to_value(u.Dalton)
        #sum flux over flavors
        array_total = sum([flux[flv].array for flv in self.flavor])
        #create a summary flux container
        flux_total = Container(array_total, flavor=self.flavor, 
                               time=flux.time, energy=flux.energy,
                               integrable_axes=flux._integrable_axes)
        rate = self.xsec*flux_total*self.weight*Ntargets
        return rate

class Detector:
    """A detector configuration for the rate calculation. """
    
    def __init__(self, name:str, mass: u.Quantity, channels: Dict[str,DetectionChannel]):
        """
        Parameters
        ----------
        name: str
            Detector name
        mass: Quantity[mass]
            Detector mass
        channels: Dict[str,DetectionChannel]
            Dictionary of detection channels in the format {name:channel}
    
        Note
        ----
        These parameters and detection channels can be modified later, before calling :meth:`run`
        """
        self.name = name
        self.mass = mass
        self.channels = channels

    def __repr__(self):
        return f'Detector(name="{self.name}", mass={self.mass}, channels={list(self.channels)})'
    
    def run(self, flux:Container, detector_effects:bool=True)->Dict[str, Container]:
        """Calculate the interaction rates for all channels in the detector.

        Parameters
        ----------
        flux:Container
            The incoming neutrino flux (or fluence).
        detector_effects:bool
            If `True` (default) apply the smearing and efficiency for all channels. Otherwise just calculate the interaction rates vs neutrino energy.

        Returns
        -------
        Dict[str,Container]
            Event rate for each detection channel in as dictionary {channel name: event rate}
        """
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
        """Read the detector configuration from the SNOwGLoBES

        Parameters
        ----------
        name:str
            Detector name (see :attr:`detectors` for options)
        material:str or None
            Detector material (see :attr:`materials` for options)
            If `None` (default) try to guess material from detector name

        Returns
        -------
        Detector
            an object with the detector configuration.
        """
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
            channel = DetectionChannel(
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