from typing import Union,Optional,Set
from snewpy.neutrino import Flavor
from astropy import units as u

import numpy as np

from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d
from enum import IntEnum

from copy import copy
from functools import wraps

class Axes(IntEnum):
    """Number of the array dimension for each axis""" 
    flavor=0
    time=1
    energy=2
    @classmethod
    def get(cls, value:Union['Axes',int,str])->'Axes':
        "convert string,int or Axes value to Axes"
        if isinstance(value,str):
            return cls[value]
        else:
            return cls(value)

class _ContainerBase:
    unit = None
    def __init__(self, 
                 data: u.Quantity,
                 flavor: np.array,
                 time: u.Quantity[u.s], 
                 energy: u.Quantity[u.MeV],
                 *,
                 integrable_axes: Optional[Set[Axes]] = None
    ):
        if self.unit is not None:
            #try to convert to the unit
            data = data.to(self.unit)
        self.array = data
        self.flavor = np.array(flavor, subok=True)
        self.time = time
        self.energy = energy
        
        if integrable_axes is not None:
            #store which axes can be integrated
            self._integrable_axes = set(integrable_axes)
        else:
            #guess which axes can be integrated
            self._integrable_axes = {a for a in Axes if(self._axshape[a]==self.array.shape[a])}
            self._integrable_axes.discard(Axes.flavor)
            
    @property
    def _sumable_axes(self):
        return set(Axes).difference(self._integrable_axes)
    @property
    def axes(self):
        return self.flavor, self.time, self.energy

    @property
    def _axshape(self):
        return tuple(len(a) for a in self.axes)

    @property
    def shape(self):
        return self.array.shape
        
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
        newaxes = [ax.__getitem__(arg) for arg, ax in zip(args, self.axes)]
        return self.__class__(array, *newaxes)

    def __repr__(self) -> str:
        """print information about the container"""
        s = [f"{len(values)} {label.name}({values.min()};{values.max()})"
            for label, values in zip(Axes,self.axes)]
        return f"{self.__class__.__name__} {self.array.shape} [{self.array.unit}]: <{' x '.join(s)}>"
    
    def sum(self, axis: Union[Axes,str])->'Container':
        """sum along given axis, producing a reduced array"""
        axis = Axes.get(axis)
        if axis not in self._sumable_axes:
            raise ValueError(f'Cannot sum over {axis.name}! Valid axes are {self._sumable_axes}')
        array = np.sum(self.array, axis = axis, keepdims=True)
        axes = list(self.axes)
        axes[axis] = axes[axis].take([0,-1])
        return Container(array,*axes, integrable_axes = self._integrable_axes.difference({axis}))

    def integrate(self, axis:Union[Axes,str], limits:np.ndarray=None)->'Container':
        """integrate along given axis, producing a reduced array"""
        axis = Axes.get(axis)
        if not axis in self._integrable_axes:
            raise ValueError(f'Cannot integrate over {axis.name}! Valid axes are {self._integrable_axes}')
        #set the limits
        ax = self.axes[axis]
        xmin, xmax = ax.min(), ax.max()
        if limits is None:
            limits = u.Quantity([xmin, xmax])
        limits = limits.to(ax.unit)
        limits = limits.clip(xmin, xmax)
        #compute the integral
        yc = cumulative_trapezoid(self.array, x=ax, axis=axis, initial=0)
        _integral = interp1d(x=ax, y=yc, fill_value=0, axis=axis, bounds_error=False)
        array = np.diff(_integral(limits),axis=axis) << (self.array.unit*ax.unit)
        axes = list(self.axes)
        axes[axis] = limits
        #choose the proper class
        return Container(array, *axes, integrable_axes=self._integrable_axes.difference({axis}))
        return result
        
    def integrate_or_sum(self, axis:Union[Axes,str])->'Container':
        if self.can_integrate(axis):
            return self.integrate(axis)
        else:
            return self.sum(axis)
            
    def can_integrate(self, axis):
         return Axes.get(axis) in self._integrable_axes
    def can_sum(self, axis):
         return Axes.get(axis) not in self._integrable_axes
    
    def __rmul__(self, factor):
        "multiply array by givem factor or matrix"
        return self.__mul__(self, factor)

    def __mul__(self, factor) -> 'Container':
        "multiply array by givem factor or matrix"
        #if not (np.isscalar(factor)):
        #    raise ValueError("Factor should be a scalar value")
        array = self.array*factor
        axes = list(self.axes)
        return Container(array, *axes)

    def save(self, fname:str):
        """Save container data to a NPZ file"""
        def _save_quantity(name):
            values = self.__dict__[name]
            try:
                array = values.to_value()
                unit = values.unit.to_string()
                return {name:array, f'_{name}_unit':unit}
            except:
                return {name:values}
        data_dict = {}
        for name in ['array','time','energy','flavor']:
            data_dict.update(_save_quantity(name))
        np.savez(fname,
                 _class_name=self.__class__.__name__, 
                 **data_dict,
                 _integrable_axes=np.array([int(a) for a in self._integrable_axes])
                )
    
    @classmethod
    def load(cls, fname:str)->'Container':
        """Load container from a given file"""
        with np.load(fname) as f:
            def _load_quantity(name):
                array = f[name]
                try:
                    unit = str(f[f'_{name}_unit'])
                    return array<<u.Unit(unit)
                except KeyError:
                    return array
                    
            class_name = str(f['_class_name'])
            array = _load_quantity('array')
            cls = Container[array.unit, class_name]
            return cls(data=array,
                       **{name:_load_quantity(name) for name in ['time','energy','flavor']},
                       integrable_axes=f['_integrable_axes'])
            
    def __eq__(self, other:'Container')->bool:
        "Check if two Containers are equal"
        result = self.__class__==other.__class__ and \
                 self.unit == other.unit and \
                 np.allclose(self.array, other.array) and \
                 all([np.allclose(self.axes[ax], other.axes[ax]) for ax in Axes])
        return result

class Container(_ContainerBase):
    """Choose appropriate container class for the given unit"""
    #a dictionary holding classes for each unit
    _unit_classes = {}

    @wraps(_ContainerBase.__init__)
    def __new__(cls, data,*args, **kwargs):
        if not cls.unit:
            data = u.Quantity(data)
            cls = cls[data.unit]
        return super().__new__(cls)

    def __class_getitem__(cls, args):
        try:
            unit,name = args
        except:
            unit,name = args, None
        #convert units to astropy.Unit
        unit = u.Unit(unit)
        if unit not in cls._unit_classes:
            #create subclass of this type with given unit
            name = name or f'{cls.__name__}[{unit}]'
            cls._unit_classes[unit] = type(name,(cls,),{'unit':unit})
        return cls._unit_classes[unit]
        

#some standard container classes that can be used for 
Flux = Container['1/(MeV*s*cm**2)', "d2FdEdT"]
Fluence = Container[Flux.unit*u.s, "dFdE"]
Spectrum= Container[Flux.unit*u.MeV, "dFdT"]
IntegralFlux= Container[Flux.unit*u.s*u.MeV, "dF"]

DifferentialEventRate = Container['1/(MeV*s)', "d2NdEdT"]
EventRate = Container['1/s', "dNdT"]
EventSpectrum = Container['1/MeV', "dNdE"]
EventNumber = Container['', "N"]
