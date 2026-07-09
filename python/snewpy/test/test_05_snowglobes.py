import pytest
from pathlib import Path
from snewpy.test._rate_crosscheck_table import rate_table
from snewpy import snowglobes, model_path

from snewpy.models.ccsn import *
import astropy.units as u

pytestmark=pytest.mark.snowglobes

#get available model parameters from table
param_values = list(rate_table.keys())

#get available detectors from table
detectors = list(list(rate_table.values())[0].keys())

#make sure the model files are loaded
model = Bollig_2016
for params in model.get_param_combinations():
    model(**params)
    
def fluence_calculation(model_name,model_mass,transform):
    #generating fluence file
    model = snowglobes.get_model_class(model_name)(model_mass<<u.Msun)    
    return snowglobes.generate(model, transform,d=10*u.kpc)

def rates_calculation(fluence):
    table = snowglobes.calculate(None,fluence,detector=detectors)
    result = {}
    for det in table:
        result[det] = 0
        for chan in table[det]:
            result[det] += table[det][chan].integrate_or_sum('energy').array.squeeze().value
    return result

@pytest.mark.parametrize('model_parameters',param_values)
def test_total_rate_equals_table_value(model_parameters):
    fluence = fluence_calculation(*model_parameters)
    calculated_rates  = rates_calculation(fluence)
    for detector in detectors:
        expected = pytest.approx(rate_table[model_parameters][detector], rel=0.01)
        assert calculated_rates[detector] == expected, f"Crosscheck failed for {detector}"
