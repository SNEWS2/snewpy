"""
The module ``snewpy.rate_calculator`` contains a Python interface for calculating event rates using data from SNOwGLoBES

"""
import numpy as np
from snewpy.snowglobes_interface import SnowglobesData, guess_material
from snewpy.neutrino import Flavor
from snewpy.flux import Container
from astropy import units as u
from warnings import warn

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
    
#--------------------------------------
class RateCalculator(SnowglobesData):
    def __init__(self, base_dir=''):
        super().__init__(base_dir=base_dir)
        
    def calc_rate(self, channel, TargetMass, flux):
        """calculate interaction rate"""
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
    
    def run(self, flux, detector:str, material:str=None, detector_effects:bool = True):
        TargetMass = float(self.detectors[detector].tgt_mass)<<(1e3*u.tonne)
        material = material or guess_material(detector)
        binning = self.binning[material]
        energies_t = binning['e_true']<<u.GeV
        energies_s = binning['e_smear']<<u.GeV
        smearing_shape = (len(energies_t)-1, len(energies_s)-1)
        result = {}
        for channel in self.channels[material].itertuples():
            rate = self.calc_rate(channel, TargetMass, flux)
            #apply channel weight 
            rate = rate*channel.weight
            
            if detector_effects:
                if rate.can_integrate('energy'):
                    #integrate over given energy bins
                    rateI = rate.integrate('energy', energies_t)
                else:
                    if not np.all(rate.energy.shape==energies_s.shape) \
                    or not np.allclose(rate.energy,energies_s):
                        raise ValueError(f'Fluence energy values should be equal to smearing matrix binning. Check binning[{material}]["e_true"]!')
                    rateI = rate
                try:
                    smear = self.smearings[detector][channel.name]
                    effic = self.efficiencies[detector][channel.name]
                except KeyError:
                    warn(f'Detector effects not found for detector={detector}, channel={channel.name}. Using unsmeared and with 100% efficiency')
                    smear = np.eye(*smearing_shape)
                    effic = np.ones(len(energies_s)-1)
                #apply smearing
                rateS_array = np.dot(rateI.array, smear.T) * effic
                rateS = Container(rateS_array, rateI.flavor, rateI.time, energies_s)
                rateI = rateS
                
                #result[(channel.name,'unsmeared','unweighted')] = rateI
                result[channel.name] = rateS
            else:
                result[channel.name] = rate
        return result

    def run_simple(self, flux, detector:str, material:str=None, detector_effects:bool = True):
        """This function tries to reproduce the calculations in generate_fluence/simulate
        Main difference: use bin centers to evaluate the cross-sections
        """
        TargetMass = float(self.detectors[detector].tgt_mass)<<(1e3*u.tonne)
        material = material or guess_material(detector)
        binning = self.binning[material]
        energies_t = binning['e_true']<<u.GeV
        energies_s = binning['e_smear']<<u.GeV
        binwidth = np.diff(energies_t).to(flux.energy.unit)
        
        smearing_shape = (len(energies_t)-1, len(energies_s)-1)
        #check that energy is the same
        assert np.allclose(flux.energy,center(energies_t))
        
        result = {}
        for channel in self.channels[material].itertuples():
            rate = self.calc_rate(channel, TargetMass, flux)
            #apply channel weight 
            rate = rate*channel.weight
            #"integrate" over energy bins
            rateI = rate*binwidth
            if detector_effects:
                smear = self.smearings[detector].get(channel.name, 
                                                     np.eye(*smearing_shape)
                                                    )
                effic = self.efficiencies[detector].get(channel.name, 
                                                        np.ones(len(energies_s)-1)
                                                       )
                #apply smearing
                rateS_array = np.dot(rateI.array, smear.T) * effic
                rateS = Container(rateS_array, rateI.flavor, rateI.time, energies_s)
                rateI = rateS
                
                #result[(channel.name,'unsmeared','unweighted')] = rateI
                result[channel.name] = rateI
            else:
                result[channel.name] = rate
        return result