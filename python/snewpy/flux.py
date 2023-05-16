from snewpy.neutrino import Flavor
from astropy import units as u

import numpy as np

from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d
from enum import Enum

class _Container:
    class Axes(Enum):
        flavor=0
        time=1
        energy=2
        
    def __init__(self, data: u.Quantity,
                       flavor: np.array,
                       time: u.Quantity[u.s], 
                       energy: u.Quantity[u.MeV],
                ):
        self.array = data.to(self.unit)
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
        array = np.sum(self.array, axis = axis.value, keepdims=True)
    
    def integrate(self, axis:Axes, limits:np.ndarray = None):
        #set the limits
        ax = self._axes[axis.value]
        xmin, xmax = ax.min(), ax.max()
        if limits is None:
            limits = u.Quantity([xmin, xmax])

        limits = limits.to(ax.unit)
        limits = limits.clip(xmin, xmax)

        #compute the integral
        x  = self._axes[axis.value]
        yc = cumulative_trapezoid(self.array, x=x, axis=axis.value, initial=0)
        _integral = interp1d(x=x, y=yc, fill_value=0, axis=axis.value, bounds_error=False)

        array = np.diff(_integral(limits),axis=axis.value) << (self.array.unit*ax.unit)
        axes = list(self._axes)
        axes[axis.value] = limits

        #choose the proper class
        return _choose_class(array.unit)(array, *axes)
class Events (_Container):
    unit = u.dimensionless_unscaled/u.cm**2        
class Flux (_Container):
    unit = Events.unit/u.s/u.MeV
class Fluence (_Container):
    unit = Events.unit/u.MeV
class Rate (_Container):
    unit = Events.unit/u.s

def _choose_class(unit):
    for cls in Flux,Fluence,Rate,Events:
        if cls.unit.is_equivalent(unit):
            return cls
    raise RuntimeError(f"Class not found for unit {unit}")
