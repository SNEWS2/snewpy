import numpy as np
import astropy.units as u
import pytest
from snewpy.models import ccsn
from snewpy.rate_calculator import RateCalculator

from snewpy import flavor_transformation as ft
pytestmark=pytest.mark.snowglobes

#numbers from SNEWSv2 paper - for approximate checks 
crosscheck_table=[  ('icecube',          320e3),
                    ('wc100kt30prct',4044.8/0.32),
                    ('ar40kt',              2700),
                    ('novaND',                40),
                    ('novaFD',              1900),
                    ('halo1',                  4),
                    ('halo2',                 53)]

#numbers from the SNEWPY SNOwGLoBES interface - for approximate checks 
#of simple rate computation (no smearing)
crosscheck_table_unsmeared=[  ('icecube',          7.2352e6),
                    ('wc100kt30prct',4486.9/0.32),
                    ('ar40kt',              2858.7),
                    ('novaND',                41.7231),
                    ('novaFD',              1947.0768),
                    ('halo1',                  11.5709),
                    ('halo2',                 146.4671)]

@pytest.fixture(autouse=True)
def snmodel():
    return ccsn.Bollig_2016(progenitor_mass=11.2<<u.Msun)

@pytest.fixture(autouse=True)
def rc():
    return RateCalculator()
    
@pytest.fixture
def fluence(snmodel):
    times    = snmodel.get_time()
    energies = np.linspace(0.5,100,201)<<u.MeV
    #get the flux from the model
    flux = snmodel.get_flux(t=times, E=energies, distance=10<<u.kpc, flavor_xform=ft.AdiabaticMSW())
    return flux.integrate('time')

@pytest.mark.parametrize('detector, expected_total',crosscheck_table_unsmeared)
def test_rate_unsmeared(rc, fluence, detector, expected_total):    
    #calculate the event rates
    rates = rc.run(fluence, detector, detector_effects=False)
    #calculate the total IBD event rate 
    N_total = sum([r.integrate_or_sum('energy').array.value for r in rates.values()])
    print(N_total)
    #check the final value
    assert N_total == pytest.approx(expected_total, 0.01)
    
@pytest.mark.parametrize('detector, expected_total',crosscheck_table)
def test_rate_smeared(rc, fluence, detector, expected_total):
    #calculate the event rates
    rates = rc.run(fluence, detector, detector_effects=True)
    #calculate the total IBD event rate 
    N_total = sum([r.integrate_or_sum('energy').array.value for r in rates.values()])
    #check the final value
    assert N_total == pytest.approx(expected_total, 0.1)
