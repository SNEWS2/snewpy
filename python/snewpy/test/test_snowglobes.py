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

@pytest.fixture(autouse=True)
def sng():
    return SNOwGLoBES()

@pytest.mark.parametrize('material, detector, expected_total',det_mats)
def test_snowglobes(sng,detector, material, expected_total):
    flux = './Bollig_2016_s11.2c_AdiabaticMSW_NMO.dat'
    result = {}
    material = get_material(detector)
    data = sng.run(flux,detector,material)
    total = data.weighted.smeared.sum().sum()
    assert total == pytest.approx(expected_total, 0.1)

from snewpy.snowglobes import simulate,collate, generate_time_series
from pathlib import Path

@pytest.mark.parametrize('model',['Bollig_2016/s11.2c', 'Bollig_2016/s27.0c'])
@pytest.mark.parametrize('transformation',['AdiabaticMSW_NMO','AdiabaticMSW_IMO'])
def test_generate_time_series(sng,model,transformation):
    model_path = Path('./models')/model
    model_type = model_path.parent.stem
    output_name = f'{model.replace("/","_")}_{transformation}_times_'
    generate_time_series(model_path,model_type, transformation, d=10, output_filename=output_name)

def test_simulate(benchmark):
    tarball_name='./models/Bollig_2016/Bollig_2016_s27.0c_AdiabaticMSW_IMO_times_kpc.tar.bz2'
    r = benchmark(simulate, None,tarball_name,'icecube')

def test_collate():
    tarball_name='./models/Bollig_2016/Bollig_2016_s27.0c_AdiabaticMSW_IMO_times_kpc.tar.bz2'
    collate(None, tarball_name,'icecube')
