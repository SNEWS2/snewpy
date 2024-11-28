import os
from tempfile import TemporaryDirectory

from hypothesis import given, assume, note
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays, array_shapes

from snewpy.flux import Container,Axes
from snewpy.neutrino import Flavor
import numpy as np
import astropy.units as u
import pytest

#define strategies to generate float values
def float_values(shape=array_shapes(max_dims=1),**kwargs):
    return arrays(float, shape=shape, 
                  elements={'allow_nan':False, 
                            'allow_infinity':False,
                            'min_value':-1e100,
                            'max_value':1e100
                           },
                  **kwargs)

def sorted_floats(shape=array_shapes(max_dims=1)):
    return float_values(shape, unique=True).map(np.sort)

def to_unit(unit):
    def _f(values):
        return values<<u.Unit(unit)

@st.composite
def quantities(draw, values, unit):
    return draw(values)<<unit

#how to generate axes values
flavors = st.lists(elements=st.sampled_from(Flavor),
                   unique=True,min_size=1,max_size=len(Flavor)).map(sorted)
units = st.sampled_from(['1/MeV','1/(MeV*s)','1/(MeV*s*m**2)']).map(u.Unit)
energies = quantities(sorted_floats(shape=array_shapes(max_dims=1)),u.MeV)
times    = quantities(sorted_floats(shape=array_shapes(max_dims=1)),u.s)
integrable_axess = st.lists(st.sampled_from(list(Axes)[1:]))

@st.composite
def random_flux_containers(draw, times=times, energies=energies, units=units, axes=st.just(None)):
    "Strategy to generate a random flux container"
    T,E = draw(times),draw(energies)
    flavor = draw(flavors)
    unit = draw(units)
    
    shape=[len(flavor),len(T),len(E)]
    array_data = quantities(float_values(shape),unit)
    array = draw(array_data)
    return Container[unit](array,flavor,T,E, integrable_axes=draw(axes))

#remove infinities and NaNs for the floats generation
st.register_type_strategy(float, st.floats(allow_nan=False, allow_infinity=False))

#starting test sets
#--------------------------------------------------------------
@given(flavor=flavors, time=times, energy=energies, unit=units)
def test_construction_succesfull(flavor, energy, time, unit):
    data = np.ones([len(flavor),len(time), len(energy)])<<unit
    assert Container[unit](data, flavor, time, energy)

@given(flavor=flavors, time=times, energy=energies, unit=units)
def test_construction_with_wrong_units_raises_ConversionError(flavor, energy, time, unit):
    data = np.ones([len(flavor),len(time), len(energy)])<<unit
    with pytest.raises(u.UnitConversionError):
        Container[unit*u.kg](data, flavor, time, energy)
    with pytest.raises(u.UnitConversionError):
        Container[unit](data*u.kg, flavor, time, energy)

@given(flavor=flavors, time=times, energy=energies, unit=units)
def test_construction_with_wrong_dimensions_raises_ValueError(flavor, energy, time, unit):
    data = np.ones([len(flavor)+1,len(time), len(energy)])<<unit
    with pytest.raises(ValueError):
        Container[unit](data, flavor, time, energy)

@given(f=random_flux_containers())
def test_summation_over_flavor(f:Container):
    fS = f.sum('flavor')
    #check the result array shapes
    assert fS.shape[Axes.flavor] == 1
    assert fS.shape[Axes.energy] == f.shape[Axes.energy]
    assert fS.shape[Axes.time]   == f.shape[Axes.time]
    #check the resulting array values
    assert np.allclose(fS.array, f.array.sum(axis=Axes.flavor))

axess_not_flavor = st.sampled_from([Axes.time,Axes.energy])

@given(f=random_flux_containers(), axis=axess_not_flavor
      )
def test_summation_over_wrong_axis_raises_ValueError(f:Container, axis:str):
    with pytest.raises(ValueError):
        fS = f.sum(axis)
    
@given(f=random_flux_containers())
def test_integration_over_flavor_raises_ValueError(f:Container):
    with pytest.raises(ValueError):
        fS = f.integrate('flavor')
    
@given(f=random_flux_containers(), axis=axess_not_flavor)
def test_integration_over_axes(f:Container, axis):
    #we don't want to integrate over single value
    assume(f.shape[axis]>1)
    note(f'Array is {str(f.array)}')
    fI = f.integrate(axis)
    for a in Axes:
        if a!=axis:
            assert fI.shape[a] == f.shape[a]
        else:
            assert fI.shape[a] == 1
    #check the resulting array values
    assert fI.unit == f.unit*f.axes[axis].unit

@given(f=random_flux_containers())
def test_save_and_load(f):
    with TemporaryDirectory() as tmpdir:
        fname = os.path.join(tmpdir, 'flux.npz')
        f.save(fname)
        f1 = Container.load(fname)
        assert f1 == f