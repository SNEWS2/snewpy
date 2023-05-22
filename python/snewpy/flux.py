from typing import Union
from snewpy.neutrino import Flavor
from astropy import units as u

import numpy as np

from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d
from enum import IntEnum

from copy import copy

class Axes(IntEnum):
        """Number of the array dimension for each axis""" 
        flavor=0
        time=1
        energy=2
            
class Container:
    def __init__(self, 
                 data: u.Quantity,
                 flavor: np.array,
                 time: u.Quantity[u.s], 
                 energy: u.Quantity[u.MeV],
                 *,
                 integrable_axes = None
    ):
        self.array = data
        self.flavor = flavor
        self.time = time
        self.energy = energy
        #guess which axes can be integrated
        if(integrable_axes):
            self._integrable_axes = integrable_axes
        else:
            self._integrable_axes = {a for a in Axes if(self._axshape[a]==self.array.shape[a])}
            self._integrable_axes.discard(Axes.flavor)
    @property
    def _sumable_axes(self):
        return set(Axes).difference(self._integrable_axes)
    @property
    def _axes(self):
        return self.flavor, self.time, self.energy

    @property
    def _axshape(self):
        return tuple(len(a) for a in self._axes)

    def __getitem__(self, args)->'Container':
        """Slice the flux array and produce a new Flux object"""
        try: 
            iter(args)
        except TypeError:
            args = [args]
        args = [a if isinstance(a, slice) else slice(a, a + 1) for a in args]
        #expand args to match axes
        args+=[slice(None)]*(len(Axes)-len(args))
        array = self.array.__getitem__(tuple(args))
        newaxes = [ax.__getitem__(arg) for arg, ax in zip(args, self._axes)]
        return self.__class__(array, *newaxes)

    def __repr__(self) -> str:
        """print information about the container"""
        s = [f"{len(values)} {label.name}({values.min()};{values.max()})"
            for label, values in zip(Axes,self._axes)]
        return f"{self.__class__.__name__} {self.array.shape} [{self.array.unit}]: <{' x '.join(s)}>"
    
    def sum(self, axis: Union[Axes,str])->'Container':
        """sum along given axis, producing a reduced array"""
        if isinstance(axis, str):
            axis = Axes[axis]
        if axis not in self._sumable_axes:
            raise ValueError(f'Cannot sum over {axis.name}! Valid axes are {self._sumable_axes}')
        array = np.sum(self.array, axis = axis, keepdims=True)
        axes = list(self._axes)
        axes[axis] = axes[axis].take([0,-1])
        return ContainerClass(array.unit)(array,*axes, integrable_axes = self._integrable_axes.difference({axis}))

    def integrate(self, axis:Union[Axes,str], limits:np.ndarray=None)->'Container':
        """integrate along given axis, producing a reduced array"""
        if isinstance(axis, str):
            axis = Axes[axis]
        if not axis in self._integrable_axes:
            raise ValueError(f'Cannot integrate over {axis.name}! Valid axes are {self._integrable_axes}')
        #set the limits
        ax = self._axes[axis]
        xmin, xmax = ax.min(), ax.max()
        if limits is None:
            limits = u.Quantity([xmin, xmax])
        limits = limits.to(ax.unit)
        limits = limits.clip(xmin, xmax)
        #compute the integral
        yc = cumulative_trapezoid(self.array, x=ax, axis=axis, initial=0)
        _integral = interp1d(x=ax, y=yc, fill_value=0, axis=axis, bounds_error=False)
        array = np.diff(_integral(limits),axis=axis) << (self.array.unit*ax.unit)
        axes = list(self._axes)
        axes[axis] = limits
        #choose the proper class
        return ContainerClass(array.unit)(array, *axes, integrable_axes=self._integrable_axes.difference({axis}))
        return result

    def __rmul__(self, factor):
        "multiply array by givem factor or matrix"
        return self.__mul__(self, factor)

    def __mul__(self, factor) -> 'Container':
        "multiply array by givem factor or matrix"
        #if not (np.isscalar(factor)):
        #    raise ValueError("Factor should be a scalar value")
        array = self.array*factor
        axes = list(self._axes)
        return ContainerClass(array.unit)(array, *axes)

#a dictionary holding classes for each unit
_unit_classes = {}

def ContainerClass(Unit, name=None):
    """Choose appropriate container class for the given unit"""
    if Unit in _unit_classes:
        return _unit_classes[Unit]
        
    class _cls(Container):
        unit = Unit
        def __init__(self, data: u.Quantity,
                           flavor: np.array,
                           time: u.Quantity[u.s], 
                           energy: u.Quantity[u.MeV]
                    ):
            
            super().__init__(data.to(self.unit),flavor,time,energy)
            
    if(name):
        _cls.__name__=name
    #register unit classes
    _unit_classes[Unit] = _cls
    #return the registered class
    return _unit_classes[Unit]

#some standard container classes that can be used for 
Flux = ContainerClass(u.one/u.MeV/u.s/u.cm**2, "d2FdEdT")
Fluence = ContainerClass(Flux.unit*u.s, "dFdE")
Spectrum= ContainerClass(Flux.unit*u.MeV, "dFdT")

DifferentialEventRate = ContainerClass(u.one/u.MeV/u.s, "d2NdEdT")
EventRate = ContainerClass(u.one/u.s, "dNdT")
EventSpectrum = ContainerClass(u.one/u.MeV, "dNdE")
EventNumber = ContainerClass(u.one, "N")
