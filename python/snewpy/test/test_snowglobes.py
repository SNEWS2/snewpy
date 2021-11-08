import pytest
import numpy as np
from snewpy.snowglobes import simulate,collate, generate_time_series, generate_fluence, SNOwGLoBES
from pathlib import Path
import tarfile

pytestmark=pytest.mark.snowglobes

models = ['Bollig_2016/s11.2c', 'Bollig_2016/s27.0c']
transformations = ['AdiabaticMSW_NMO','AdiabaticMSW_IMO']

@pytest.mark.parametrize('model',models)
@pytest.mark.parametrize('transformation',transformations)
def test_generate_fluence(sng,model,transformation):
    """make some input fluxes"""
    model_path = Path('./models')/model
    model_type = model_path.parent.stem
    output_name = f'fluence_{model.replace("/","_")}_{transformation}'
    fname = generate_fluence(model_path,model_type, transformation, d=10, output_filename=output_name)
    with tarfile.open(fname) as tar:
        tar.extractall(model_path.parent)

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

#numbers from SNEWSv2 paper - for approximate checks 
crosscheck_table=[  ('icecube',          320e3),
                    ('wc100kt30prct',4044.8/0.32),
                    ('ar40kt',              2700),
                    ('novaND',                40),
                    ('novaFD',              1900),
                    ('halo1',                  4),
                    ('halo2',                 53)]


@pytest.fixture(autouse=True)
def sng():
    return SNOwGLoBES()

@pytest.mark.parametrize('detector, expected_total',crosscheck_table)
def test_snowglobes_crosscheck(sng, detector, expected_total):
    flux = './models/Bollig_2016/fluence_Bollig_2016_s11.2c_AdiabaticMSW_NMO.dat'
    material = get_material(detector)
    data = sng.run(flux,detector,material)
    total = data[0].weighted.smeared.sum().sum()
    assert total == pytest.approx(expected_total, 0.1)


def process(tarball_name):
    simulate(None,tarball_name,'icecube')
    collate(None, tarball_name,'icecube')

@pytest.mark.timing    
def test_simulation_chain_benchmark(benchmark):
    tarball_name='./models/Bollig_2016/fluence_Bollig_2016_s27.0c_AdiabaticMSW_IMO.tar.bz2'
    r = benchmark(process,tarball_name)
