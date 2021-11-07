import pytest
import numpy as np
from snewpy.snowglobes_interface import SNOwGLoBES

def get_material(detector):
    if detector.startswith('wc') or detector.startswith('ice'):
        return 'water'
    elif detector.startswith('ar'):
        return 'argon'
    elif detector.startswith('nova'):
        return 'nova_soup'
    elif detector.startswith('halo'):
        return 'lead'
    elif detector.startswith('scint'):
        return 'scint'

det_mats =[ ('water','icecube',      320e3),
            ('water','wc100kt30prct',4044.8/0.32),
            ('argon','ar40kt',       2700),
            ('nova_soup','novaND',   40),
            ('nova_soup','novaFD',   1900),
            ('lead','halo1',         4),
            ('lead','halo2',         53)]

@pytest.mark.parametrize('material, detector, expected_total',det_mats)
def test_snowglobes(detector, material, expected_total):
    sng = SNOwGLoBES()
    flux = './Bollig_2016_s11.2c_AdiabaticMSW_NMO.dat'
    result = {}
    material = get_material(detector)
    data = sng.run(flux,detector,material)
    total = data.weighted.smeared.sum().sum()
    assert total == pytest.approx(expected_total, 0.1)

from snewpy.snowglobes import simulate

def test_simulate():
   r = simulate(None,'./Bollig_2016_s11.2c_AdiabaticMSW_NMO.tar.bz2','all')
   print(r)
   assert False

