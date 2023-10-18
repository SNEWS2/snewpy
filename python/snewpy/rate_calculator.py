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

#various utility methods
def center(a):
    return 0.5*(a[1:]+a[:-1])

def _get_flavor_index(channel):
    _map = {'+e':Flavor.NU_E,
            '-e':Flavor.NU_E_BAR,
            '+m':Flavor.NU_X,
            '-m':Flavor.NU_X_BAR,
            '+t':Flavor.NU_X,
            '-t':Flavor.NU_X_BAR
            }
    return _map[channel.parity+channel.flavor]

def _get_xsec_column(channel):
    _map = {'+e':1, '+m':2, '+t':3, '-e':4, '-m':5, '-t':6 }
    return _map[channel.parity+channel.flavor]
    
def _load_xsec(self, channel, energies):
    """Load cross-section for a given channel, interpolated in the energies"""
    xsec = np.loadtxt(self.base_dir/f"xscns/xs_{channel.name}.dat")
    # Cross-section in 10^-38 cm^2
    E = energies.to_value('GeV')
    xp = xsec[:,0]
    yp = xsec[:, _get_xsec_column(channel)]
    xsecs = np.interp(np.log(E)/np.log(10), xp, yp, left=0, right=0)*E*1e-38 <<u.cm**2
    return xsecs

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
        
    def _calc_rate(self, channel, TargetMass, flux):
        """calculate interaction rate for given channel"""
        flavor = _get_flavor_index(channel)
        #check if we are provided with energy bins:
        if not flux.can_integrate('energy'):
            #we have bins, let's use central values for sampling
            Evalues = center(flux.energy)
        else:
            #we have sample points, let's just use them
            Evalues = flux.energy
        xsec = _load_xsec(self, channel, Evalues)
        Ntargets = TargetMass.to_value(u.Dalton)
        rate = flux[flavor]*xsec*Ntargets
        return rate
    
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
        TargetMass = float(self.detectors[detector].tgt_mass)<<(1e3*u.tonne)
        material = material or guess_material(detector)
        binning = self.binning[material]
        energies_t = binning['e_true']<<u.GeV
        energies_s = binning['e_smear']<<u.GeV
        smearing_shape = (len(energies_t), len(energies_s))
        result = {}
        for channel in self.channels[material].itertuples():
            rate = self._calc_rate(channel, TargetMass, flux)
            #apply channel weight 
            rate = rate*channel.weight
            #check the energy binning if needed
            energy_bin_centers = center(rate.energy)
            if detector_effects and not rate.can_integrate('energy'):
                    if not np.allclose(energy_bin_centers<<u.GeV,energies_s<<u.GeV):
                        raise ValueError(f'Fluence energy values should be equal to smearing matrix binning. Check binning["{material}"]["e_true"]!')
            
            if detector_effects:
                if rate.can_integrate('energy'):
                    #integrate over given energy bins
                    energy_bins = _bin_edges_from_centers(energies_t)
                    rateI = rate.integrate('energy', energy_bins)
                else:
                    rateI = rate
                try:
                    smear = self.smearings[detector][channel.name]
                except KeyError:
                    warn(f'Smearing not found for detector={detector}, channel={channel.name}. Using unsmeared spectrum')
                    smear = np.eye(*smearing_shape)
                try:
                    effic = self.efficiencies[detector][channel.name]
                except KeyError:
                    warn(f'Efficiency not found for detector={detector}, channel={channel.name}. Using 100% efficiency')
                    effic = np.ones(len(energies_s))
                #apply smearing and efficiency
                rateS_array = np.dot(rateI.array, smear.T) * effic
                rateS = Container(rateS_array, 
                                  flavor=rateI.flavor, 
                                  time=rateI.time, 
                                  energy=_bin_edges_from_centers(energies_s)<<u.MeV,
                                  integrable_axes=rateI._integrable_axes)

                result[channel.name] = rateS
            else:
                result[channel.name] = rate
        return result
