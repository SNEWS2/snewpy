import pytest
import numpy as np
from snewpy.snowglobes import simulate,collate, generate_time_series, generate_fluence, SimpleRate
from pathlib import Path
import tarfile

pytestmark=pytest.mark.snowglobes

models = ['Bollig_2016/s11.2c', 'Bollig_2016/s27.0c']
transformations = ['AdiabaticMSW_NMO','AdiabaticMSW_IMO']

@pytest.mark.parametrize('model',models)
@pytest.mark.parametrize('transformation',transformations)
def test_generate_fluence(splr, model, transformation):
    """make some input fluxes"""
    model_path = Path('./models')/model
    model_type = model_path.parent.stem
    output_name = f'fluence_{model.replace("/","_")}_{transformation}'
    fname = generate_fluence(model_path,model_type, transformation, d=10, output_filename=output_name)
    with tarfile.open(fname) as tar:
        tar.extractall(model_path.parent)

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
def splr():
    return SimpleRate()

@pytest.mark.parametrize('detector, expected_total_unsmeared',crosscheck_table_unsmeared)
def test_simplerate_crosscheck(splr, detector, expected_total_unsmeared):
    flux = './models/Bollig_2016/fluence_Bollig_2016_s11.2c_AdiabaticMSW_NMO.dat'
    data = splr.run(flux,detector)
    total = data[0].weighted.unsmeared.sum().sum()
    assert total == pytest.approx(expected_total_unsmeared, 0.01)

@pytest.mark.parametrize('detector, expected_total',crosscheck_table)
def test_simplerate_smear_crosscheck(splr, detector, expected_total):
    flux = './models/Bollig_2016/fluence_Bollig_2016_s11.2c_AdiabaticMSW_NMO.dat'
    data = splr.run(flux,detector)
    total = data[0].weighted.smeared.sum().sum()
    assert total == pytest.approx(expected_total, 0.1)
