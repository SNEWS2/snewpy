import astropy.units as u
import numpy as np
from snewpy.models import pisn
from snewpy.neutrino import MassHierarchy, MixingParameters
from snewpy.flavor_transformation import AdiabaticMSW
from snewpy.rate_calculator import RateCalculator
import pytest

pytestmark=pytest.mark.snowglobes

rc = RateCalculator()

distance = 1000*u.pc

@pytest.mark.parametrize('model_class,model_params',[
    (pisn.Wright_2017, {'progenitor_mass': 150*u.Msun, 'eos': 'SFHo'}),
    (pisn.Wright_2017, {'progenitor_mass': 150*u.Msun, 'eos': 'Helm'}),    
    (pisn.Wright_2017, {'progenitor_mass': 250*u.Msun, 'eos': 'SFHo'}),
    (pisn.Wright_2017, {'progenitor_mass': 250*u.Msun, 'eos': 'Helm'})    
])
@pytest.mark.parametrize('transformation',[AdiabaticMSW(MixingParameters(mh)) for mh in MassHierarchy])

def test_pisn_rate(model_class, model_params, transformation):
    model = model_class(**model_params)
    times    = model.get_time()
    energies = np.linspace(0,40,201)<<u.MeV
    fluence = model.get_flux(times, energies, distance=distance, flavor_xform=transformation).integrate_or_sum('time')
    dNdE = rc.run(fluence, detector='wc100kt30prct', detector_effects=False)['ibd']
    #the factor of 0.5 is because the mass of SuperK is 50kt, not 100kt
    Nevents = 0.5 * dNdE.integrate_or_sum('energy').array.squeeze()
    assert 0.1<Nevents<10
