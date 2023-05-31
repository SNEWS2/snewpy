from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays, array_shapes

from snewpy.flux import Container,Axes
from snewpy.neutrino import Flavor
import numpy as np
import astropy.units as u
import pytest
#define strategies to generate values
def float_values(shape=array_shapes(max_dims=1)):
    return arrays(float, shape=shape, elements={'allow_nan':False, 'allow_infinity':False})

def sorted_floats(shape=array_shapes(max_dims=1)):
    return float_values(shape).map(np.sort)

def to_unit(unit):
    def _f(values):
        return values<<u.Unit(unit)

@st.composite
def quantities(draw, values, unit):
    return draw(values)<<unit

#how to generate axes values
flavors = st.lists(elements=st.sampled_from(Flavor),
                   unique=True,min_size=1,max_size=len(Flavor)).map(sorted)
units = st.sampled_from(['1/MeV','1/(MeV*s)','1/(MeV*s*cm**2)']).map(u.Unit)
energies = quantities(sorted_floats(shape=array_shapes(max_dims=1)),u.MeV)
times    = quantities(sorted_floats(shape=array_shapes(max_dims=1)),u.s)
integrable_axess = st.lists(st.sampled_from(list(Axes)[1:]))

@st.composite
def random_flux_containers(draw, times=times, energies=energies, units=units, axes=integrable_axess):
    "Strategy to generate a random flux container"
    E = draw(energies)
    T = draw(times)
    flavor = draw(flavors)
    unit = draw(units)
    
    shape=[len(flavor),len(T),len(E)]
    array_data = quantities(float_values(shape),unit)
    array = draw(array_data)
    return Container[unit](array,flavor,T,E)

@given(flavor=flavors, time=times, energy=energies, unit=units)
def test_construction_succesfull(flavor, energy, time, unit):
    data = np.ones([len(flavor),len(time), len(energy)])<<unit
    assert Container[unit](data, flavor, time, energy)

@given(flavor=flavors, time=times, energy=energies, unit=units)
def test_construction_with_wrong_units_raises_ValueError(flavor, energy, time, unit):
    data = np.ones([len(flavor),len(time), len(energy)])<<unit
    with pytest.raises(u.UnitConversionError):
        Container[unit*u.kg](data, flavor, time, energy)
    with pytest.raises(u.UnitConversionError):
        Container[unit](data*u.kg, flavor, time, energy)
