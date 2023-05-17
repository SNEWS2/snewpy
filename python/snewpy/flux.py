from snewpy.neutrino import Flavor
from astropy import units as u

import numpy as np

from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d
from enum import IntEnum

from copy import copy

class Container:
    class Axes(IntEnum):
        flavor=0
        time=1
        energy=2
        
    def __init__(self, data: u.Quantity,
                       flavor: np.array,
                       time: u.Quantity[u.s], 
                       energy: u.Quantity[u.MeV],
                ):
        self.array = data
        self.flavor = flavor
        self.time = time
        self.energy = energy
        
        #if(data.shape != self._axshape):
        #    print(f'{data.shape}!={self._axshape}')
        
    @property
    def _axes(self):
        return self.flavor, self.time, self.energy
    @property
    def _axshape(self):
        return tuple(len(a) for a in self._axes)
    
    def __getitem__(self, args):
        """Slice the flux array and produce a new Flux object.
        Parameters
        ----------
        args: slice | List[slice]
            Slices indices for each dimension to slice on
        Returns
        -------
        Flux
            New Flux object with the array and all the dimensions are slices of initial
        """
        try: 
            iter(args)
        except TypeError:
            args = [args]
        args = [a if isinstance(a, slice) else slice(a, a + 1) for a in args]
        array = self.array.__getitem__(tuple(args))
        newaxes = [ax.__getitem__(arg) for arg, ax in zip(args, self._axes)]
        return self.__class__(array, *newaxes)
    
    def __repr__(self):
        s = [
            f"{len(values)} {label.name}({values.min()};{values.max()})"
            for label, values in zip(self.Axes,self._axes)
        ]
        return f"{self.__class__.__name__} {self.array.shape} [{self.array.unit}]: <{' x '.join(s)}>"
    
    def sum(self, axis: Axes):
        array = np.sum(self.array, axis = axis, keepdims=True)
        axes = list(self._axes)
        axes[axis] = axes[axis].take([0,-1])
        return ContainerClass(array.unit)(array,*axes)
    
    def integrate(self, axis:Axes, limits:np.ndarray = None):
        #set the limits
        ax = self._axes[axis.value]
        xmin, xmax = ax.min(), ax.max()
        if limits is None:
            limits = u.Quantity([xmin, xmax])

        limits = limits.to(ax.unit)
        limits = limits.clip(xmin, xmax)

        #compute the integral
        x  = self._axes[axis]
        yc = cumulative_trapezoid(self.array, x=x, axis=axis, initial=0)
        _integral = interp1d(x=x, y=yc, fill_value=0, axis=axis, bounds_error=False)

        array = np.diff(_integral(limits),axis=axis) << (self.array.unit*ax.unit)
        axes = list(self._axes)
        axes[axis] = limits

        #choose the proper class
        return _choose_class(array.unit)(array, *axes)

    def __rmul__(self, factor):
        return self.__mul__(self, factor)
        
    def __mul__(self, factor):
        #if not (np.isscalar(factor)):
        #    raise ValueError("Factor should be a scalar value")
        array = self.array*factor
        axes = list(self._axes)
        return _choose_class(array.unit)(array, *axes)


_unit_classes = {}

def ContainerClass(Unit, name=None):
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

per_unit_area = u.one/u.cm**2
Flux = ContainerClass(u.one/u.MeV/u.s/u.cm**2, "d2FdEdT")
Fluence = ContainerClass(Flux.unit*u.s, "dFdT")
Spectrum= ContainerClass(Flux.unit*u.MeV, "dFdE")

EventSpectrumRate = ContainerClass(u.one/u.MeV/u.s, "d2NdEdT")
EventRate = ContainerClass(u.one/u.s, "dNdT")

def _choose_class(unit):
    return _unit_classes.get(unit, Container)
