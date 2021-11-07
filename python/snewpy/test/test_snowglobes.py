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

det_mats =[ ('water','icecube',      123),
            ('water','wc100kt15prct',123),
            ('argon','ar40kt',       123),
            ('nova_soup','novaND',   123),
            ('nova_soup','novaFD',   123),
            ('lead','halo1',         123),
            ('lead','halo2',         123)]

@pytest.mark.parametrize('material, detector, expected_total',det_mats)
def test_snowglobes(detector, material, expected_total):
    sng = SNOwGLoBES()
    flux = './test_Kuroda_2020_AdiabaticMSW_NMO.dat'
    result = {}
    material = get_material(detector)
    data = sng.run(flux,detector,material)
    total = data.weighted.smeared.sum().sum()
    assert total == pytest.approx(expected_total)
