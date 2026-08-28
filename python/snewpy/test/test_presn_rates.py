import astropy.units as u
import numpy as np
from snewpy.models import presn
from snewpy.neutrino import MassHierarchy, MixingParameters
from snewpy.flavor_transformation import AdiabaticMSW
from snewpy.rate_calculator import RateCalculator
import pytest

pytestmark=pytest.mark.snowglobes

rc = RateCalculator()

distance = 200*u.pc
 #SNOwGLoBES detector for water Cerenkov
T = np.geomspace(-1*u.hour, -1*u.min,1000)
E = np.linspace(0,20,100)*u.MeV

@pytest.mark.parametrize('model_class,model_params',[
    (presn.Odrzywolek_2010, {'progenitor_mass': 15*u.Msun}),
    (presn.Kato_2017, {'progenitor_mass': 15*u.Msun}),
    (presn.Patton_2017, {'progenitor_mass': 15*u.Msun}),
    (presn.Yoshida_2016, {'progenitor_mass': 15*u.Msun}),
    (presn.Myers_2026, {'progenitor_mass': 15*u.Msun, 'eta': 0.8, 'bounds': 'core'}),
])
@pytest.mark.parametrize('transformation',[AdiabaticMSW(MixingParameters(mh)) for mh in MassHierarchy])
@pytest.mark.parametrize('detector', ["wc100kt30prct"])
def test_presn_rate(model_class, model_params, transformation, detector):
    model = model_class(**model_params)
    flux = model.get_flux(T, E, distance=distance, flavor_xform=transformation)
    rate = rc.run(flux, detector='scint20kt', detector_effects=False)['ibd']
    ibd_events = rate.integrate_or_sum('time').integrate_or_sum('energy').array.squeeze().value
    assert 10<ibd_events<1000
