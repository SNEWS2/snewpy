"""
The module :mod:`snewpy.flux` defines the class :class:`snewpy.flux.Container` - a container for the neutrino flux, fluence, event rates etc.

This object wraps a 3D array, and its dimensions: `flavor`, `time` and `energy`.

Usage
-----
The flux container will be produced by the :meth:`SupernovaModel.get_flux`, and consumed by the :meth:`RateCalculator.run`

A example of usage:

.. code-block:: python

    import astropy.units as u
    import numpy as np
    from snewpy.models import ccsn
    from snewpy.flavor_transformation import AdiabaticMSW
    
    #get the suernova model
    model = ccsn.Bollig_2016(progenitor_mass=27<<u.Msun)
    #define the sampling points
    energy = np.linspace(0,100,51)<<u.MeV
    time = np.linspace(0,10,1000)<<u.s
    
    #calculate the flux
    flux = model.get_flux(time, energy, distance=10<<u.kpc, flavor_xform=AdiabaticMSW())

    #optionally: integrate the flux over time bins to obtain the fluence
    fluence = flux.integrate('time', limits=np.linspace(0,1,21)<<u.s)
    
    #calculate the event rates using the RateCalculator
    from snewpy.rate_calculator import RateCalculator
    rc = RateCalculator()
    event_rates = rc.run(fluence, detector='icecube')
    ibd_rate = event_rates['ibd'] #will be also a Container instance
    #get the total values vs. time bin
    N_ibd_vs_T = ibd_rate.sum('energy')
    #or total values vs. energy bin
    N_ibd_vs_E = ibd_rate.sum('time')
    #or total number of interactions
    N_ibd_tot = ibd_rate.sum('time').sum('energy')
    #to retrieve the number we can access the `array` member
    print(N_ibd_tot.array.squeeze()) #312249.3525844854

Reference
---------

.. autoclass:: Container
    :inherited-members:

.. autoclass:: Axes

"""
from typing import Union
# from snewpy.neutrino import Flavor
from snewpy.flavor import FlavorScheme, FlavorMatrix
from astropy import units as u

import numpy as np

from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d
from enum import IntEnum
from functools import wraps
import snewpy.utils
#list of units which will be used as units for decomposition inside the Container
snewpy_unit_bases = [u.MeV, u.m, u.s, u.kg]

import matplotlib.pyplot as plt

class Axes(IntEnum):
    """Enum to keep the number number of the array dimension for each axis""" 
    flavor=0, #Flavor dimension
    time=1,   #Time dimension
    energy=2, #Energy dimension
    
    @classmethod
    def get(cls, value: Union['Axes', int, str])->'Axes':
        "convert string, int or Axes value to Axes"
        if isinstance(value,str):
            return cls[value]
        else:
            return cls(value)

def _derive_flavor_scheme(flavor:FlavorScheme|list[FlavorScheme])->FlavorScheme:
    """Obtain the flavor scheme if provided list of flavors"""
    if isinstance(flavor, type):#issubclass without isinstance(type) check raises TypeError
        if issubclass(flavor, FlavorScheme): 
            return flavor
    else:
        #get schemes from the data
        flavor_schemes = set(f.__class__ for f in flavor)
        if len(flavor_schemes)!=1:
            raise ValueError(f"Flavors {flavor} must be from a single flavor scheme, but are from {flavor_schemes}")
        else:
            return flavor_schemes.pop()
            
class _ContainerBase:
    """base class for internal use
    :noindex:
    """
    unit = None
    def __init__(self,
                 data: u.Quantity,
                 flavor: list[FlavorScheme],
                 time: u.Quantity[u.s],
                 energy: u.Quantity[u.MeV],
                 *,
                 integrable_axes: set[Axes] | None = None,
                 flavor_scheme: FlavorScheme | None = None
    ):
        """A container class storing the physical quantity (flux, fluence, rate...), which depends on flavor, time and energy.

        Parameters
        ----------
        data: :class:`astropy.Quantity`
            3D array of the stored quantity, must have dimensions compatible with (flavor, time, energy)
        
        flavor: list or a single value of :class:`snewpy.neutrino.Flavor`
            array of flavors (should be ``len(flavor)==data.shape[0]``
        
        time: :class:`astropy.Quantity`
            sampling points in time (then ``len(time)==data.shape[1]``) 
            or time bin edges (then ``len(time)==data.shape[1]+1``) 
    
        energy: :class:`astropy.Quantity`
            sampling points in energy (then ``len(energy)=data.shape[2]``) 
            or energy bin edges (then ``len(energy)=data.shape[2]+1``) 
    
        integrable_axes: set of :class:`Axes` or None
            List of axes which can be integrated.
            If None (default) this set will be derived from the axes shapes
            
        flavor_scheme: a subclass of :class:`snewpy.flavor.FlavorSchemes` or None
            A class which lists all the allowed flavors. 
            If None (default) this value will be retrieved from the ``flavor`` arguemnt.
        """
        if self.unit is not None:
            #try to convert to the unit
            data = data.to(self.unit)
        #convert the input values to arrays if they are scalar
        self.array = u.Quantity(data)
        self.time = u.Quantity(time, ndmin=1)
        self.energy = u.Quantity(energy, ndmin=1)
        self.flavor = np.array(flavor,ndmin=1, dtype=object)
        self.flavor_scheme = flavor_scheme
        if not flavor_scheme:
            self.flavor_scheme = _derive_flavor_scheme(flavor)
            
        Nf,Nt,Ne = len(self.flavor), len(self.time), len(self.energy)
        #list all valid shapes of the input array
        expected_shapes=[(nf,nt,ne) for nf in (Nf,Nf-1) for nt in (Nt,Nt-1) for ne in (Ne,Ne-1)]
        #treat special case if data is 1d array
        if self.array.ndim==1:
            #try to reshape the array to expected shape
            for expected_shape in expected_shapes:
                if np.prod(expected_shape)==self.array.size:
                    self.array = self.array.reshape(expected_shape)
                    break
        #validate the data array shape
        if self.array.shape not in expected_shapes:
            raise ValueError(f"Data array of shape {data.shape} is inconsistent with any valid shapes {expected_shapes}")
        
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
        if not isinstance(args,tuple):
            args = [args]
        args = list(args)
        arg_slices = [slice(None)]*len(Axes)
        if isinstance(args[0],str) or isinstance(args[0],FlavorScheme):
            args[0] = self.flavor_scheme[args[0]]
        for n,arg in enumerate(args):
            if isinstance(arg, int):
                arg = slice(arg, arg + 1)
            arg_slices[n] = arg

        array = self.array.__getitem__(tuple(arg_slices))
        newaxes = [ax.__getitem__(arg) for arg, ax in zip(arg_slices, self.axes)]
        return self.__class__(array, *newaxes, flavor_scheme=self.flavor_scheme)

    def __repr__(self) -> str:
        """print information about the container"""
        s = [f"{len(values)} {label.name}({values.min()};{values.max()})"
             if label!=Axes.flavor
             else f"{len(values)} {label.name}[{self.flavor_scheme}]({values.min()};{values.max()})"
             for label, values in zip(Axes,self.axes)
        ]
        return f"{self.__class__.__name__} {self.array.shape} [{self.array.unit}]: <{' x '.join(s)}>"
    
    def sum(self, axis: Axes | str)->'Container':
        """Sum along given axis, producing a Container with the summary quantity.
        
        Parameters
        -----------
            axis: :class:`Axes` or str
                An axis to sum over. String should be one of ``"flavor"``, ``"time"`` or ``"energy"``  (check :meth:`Container.can_sum`)
        Returns
        --------
            Container with summed value
                
        Raises
        ------
            ValueError 
                if the given axis cannot be summed over (i.e. it must be integrated instead, see :meth:`Container.integrate`).
                One can check summable axes with :attr:`Container._sumable_axes`

        Example
        -------
        The resulting data array will be 3D array, but the dimension, corresponding to `axis` parameter will be reduced to 1. 
        
        For an examplar Container ``a`` of a given shape::
        
            >>> a.shape
            (4, 10, 20)
            >>> a.sum('flavor').shape
            (1, 10, 20)
        
        The axis in the class will also be modified, keeping only the first and last points of the summation::
        
            >>> a.flavor
            array([0, 1, 2, 3])
            >>> a.sum('flavor').flavor
            array([0, 3])
            
        All the other axes and dimensions will be kept the same
        """
        axis = Axes.get(axis)
        if axis not in self._sumable_axes:
            raise ValueError(f'Cannot sum over {axis.name}! Valid axes are {self._sumable_axes}')
        array = np.sum(self.array, axis = axis, keepdims=True)
        axes = list(self.axes)
        axes[axis] = axes[axis].take([0,-1])
        return Container(array,*axes, integrable_axes = self._integrable_axes.difference({axis}))

    def integrate(self, axis: Axes | str, limits:np.ndarray=None)->'Container':
        """Integrate along given axis, producing a Container with the integral quantity.
        
        Parameters
        -----------
            axis: :class:`Axes` or str
                An axis to integrate over. String should be one of ``"time"`` or ``"energy"`` (check :meth:`Container.can_integrate`)

            limits: np.ndarray or None
                A sorted array (or `astropy.Quantity` consistent with the units of given axis) of integration limits
                If limits are None (default), then integrate over the the whole range of this axis
                Otherwise limits are treated as bin edges - the integration  happens within each bin limits
                
        Returns
        --------
            Container with integrated value

        Raises
        ------
            ValueError 
                if the given axis cannot be integrated over (i.e. it must be summed instead, see :meth:`Container.sum`).
        Example
        -------
        The resulting data array will be 3D array, but the dimension, corresponding to `axis` parameter will be reduced to ``len(limits)-1``. 
        
        For an examplar Container ``a`` of a given shape::
        
            >>> a.shape
            (4, 10, 20)
            >>> a.integrate('time').shape
            (4, 1, 20)
            >>> a.integrate('time', limits=[0, 0.5, 1]<<u.s).shape
            (4, 2, 20)

        The axis in the class will also be modified, keeping only the integration limits::
        
            >>> a.time
            [0. 1. 2. 3. 4. 5. 6. 7. 8. 9.] s
            >>> a.integrate('time').time
            [0., 9.] s
            >>> a.integrate('time', limits=[0, 0.5, 1]<<u.s).time
            [0., 0.5, 1.] s
        
        All the other axes and dimensions will be kept the same
        """
        axis = Axes.get(axis)
        if not axis in self._integrable_axes:
            raise ValueError(f'Cannot integrate over {axis.name}! Valid axes are {self._integrable_axes}')
        #set the limits
        ax = self.axes[axis]
        xmin, xmax = ax.min(), ax.max()
        if limits is None:
            limits = u.Quantity([xmin, xmax])
        limits = limits.to(ax.unit)
        #limits = limits.clip(xmin,xmax)
        #compute the integral
        yc = cumulative_trapezoid(self.array, x=ax, axis=axis, initial=0)
        #get first and last value to use as the fill values
        yc_limits = (yc.take(0,axis=axis), yc.take(-1,axis=axis)) 
        #this will make the _integral constant if it gets out of bounds,
        # i.e. effectively the flux outside of bounds is zero
        _integral = interp1d(x=ax, y=yc, fill_value=yc_limits, axis=axis, bounds_error=False)
        array = np.diff(_integral(limits),axis=axis) << (self.array.unit*ax.unit)
        axes = list(self.axes)
        axes[axis] = limits
        #choose the proper class
        return Container(array, *axes, integrable_axes=self._integrable_axes.difference({axis}))

    def integrate_or_sum(self, axis: Axes | str)->'Container':
        if self.can_integrate(axis):
            return self.integrate(axis)
        else:
            return self.sum(axis)
            
    def can_integrate(self, axis):
        "return true if can be integrated along given axis"
        return Axes.get(axis) in self._integrable_axes
    def can_sum(self, axis):
        "return true if can be summed along given axis"
        return Axes.get(axis) not in self._integrable_axes
    
    def __rmul__(self, factor):
        "multiply array by given factor or matrix"
        return self.__mul__(factor)

    def __mul__(self, factor) -> 'Container':
        "multiply array by given factor or matrix"
        if not (np.isscalar(factor) or isinstance(factor, np.ndarray)):
            return NotImplemented
        #    raise ValueError("Factor should be a scalar value")
        array = self.array*factor
        axes = list(self.axes)
        return Container(array, *axes)
        
    def save(self, fname:str)->None:
        """Save container data to a given file (using `numpy.savez`)"""
        def _save_quantity(name):
            values = self.__dict__[name]
            try:
                array = values.to_value()
                unit = values.unit.to_string()
                return {name:array, f'_{name}_unit':unit}
            except:
                return {name:values}
        data_dict = {}
        for name in ['array','time','energy']:
            data_dict.update(_save_quantity(name))
        data_dict['flavor'] = np.array(self.flavor, dtype=object)
        np.savez(fname,
                 _class_name=self.__class__.__name__, 
                 **data_dict,
                 _integrable_axes=np.array([int(a) for a in self._integrable_axes])
                )
    
    @classmethod
    def load(cls, fname:str)->'Container':
        """Load container from a given file"""
        with np.load(fname, allow_pickle=True) as f:
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
                 self.flavor_scheme==other.flavor_scheme and \
                 len(self.flavor)==len(other.flavor) and \
                 all(self.flavor==other.flavor) and \
                 all([np.allclose(self.axes[ax], other.axes[ax]) for ax in list(Axes)[1:]])
        return result

    def _is_full_flavor(self):
        return all(self.flavor==list(self.flavor_scheme))
                   
    def convert_to_flavor(self, flavor:FlavorScheme):
        if(self.flavor_scheme==flavor):
            return self
        return (self.flavor_scheme>>flavor)@self
    def __rshift__(self, flavor:FlavorScheme):
        return self.convert_to_flavor(flavor)
        
    def __rmatmul__(self, matrix:FlavorMatrix):
        """Multiply this flux by a FlavorMatrix"""
        if not self._is_full_flavor():
            raise RuntimeError(f"Cannot multiply flavor matrix object {self}, expected {len(self.flavor_scheme)} flavors")
        if matrix.flavor_in!=self.flavor_scheme:
            raise ValueError(f"Cannot multiply flavor matrix {matrix} by {self} - flavor scheme mismatch!")
        #apply the multiplication:
        #first add the missing dimensions for the matrix (if needed)
        f = self.array
        m = snewpy.utils.expand_dimensions_to(matrix.array, 
                                              ndim=f.ndim+1)
        
        #do the multiplication
        array = np.einsum('ij...,j...->i...',m,f)
        return Container(array, flavor=matrix.flavor_out, time=self.time, energy=self.energy)
                  
class Container(_ContainerBase):
    #a dictionary holding classes for each unit
    _unit_classes = {}

    @wraps(_ContainerBase.__init__)
    def __new__(cls, data,*args, **kwargs):
        data = data.decompose(snewpy_unit_bases) #simplify the units, reducing to the bases
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

    @classmethod
    def from_dict(cls, 
                  data_dict: dict[FlavorScheme,np.ndarray],
                  time: u.Quantity[u.s],
                  energy: u.Quantity[u.MeV],
                  *,
                  flavor_scheme: FlavorScheme | None = None,
                  integrable_axes: set[Axes] | None = None
                 ):
        """Create a new Container from given dictionary of type {Flavor: data}"""
        flavor = list(data_dict.keys())
        if flavor_scheme is None:
            flavor_scheme = _derive_flavor_scheme(flavor)
        array = np.stack([data_dict[flv] for flv in flavor_scheme])
        #check if we need to expand the dimensions
        if time.size==1 and array.ndim<3:
            array = np.expand_dims(array, axis=Axes['time'])
        if energy.size==1 and array.ndim<3:
            array = np.expand_dims(array, axis=Axes['energy'])
            
        return cls(array, flavor, time, energy, flavor_scheme=flavor_scheme, integrable_axes=integrable_axes)
        
    def project_to(self, axis='energy', squeeze=False):
        if axis=='energy':
            fP = self.integrate_or_sum('time')
        else:
            fP = self.integrate_or_sum('energy')
        x = fP.__dict__[axis]
        if squeeze:
            return x, fP.array.squeeze().T
        else:
            return x, fP
        
        
    def plot(flux, projection='energy', styles=None, **kwargs):
        x, fP = flux.project_to(projection, squeeze=False)

        if isinstance(styles, dict):
            styles = styles.get
        elif styles==None:
            styles = lambda flv: {'ls':'-' if flv.is_neutrino else ':',
                                  'color':f'C{flv//2:d}'}
        lines = []
        for idx,flv in zip(range(fP.array.shape[0]),fP.flavor):
            style = styles(flv)
            style.update(kwargs)
            style.setdefault('label',flv.to_tex())
            y = fP[idx].array.squeeze()
            if len(x)==len(y):
                l=plt.plot(x,y, **style)
            else:
                l=plt.stairs(y,edges=x, **style)
            lines.append(l)
            
        plt.xlabel(f'{projection}, {x.unit._repr_latex_()}')
        plt.ylabel(f'{fP.__class__.__name__}, {x.unit._repr_latex_()}')
        return lines

#some standard container classes that can be used for 
Flux = Container['1/(MeV*s*m**2)', "d2FdEdT"]
Fluence = Container[Flux.unit*u.s, "dFdE"]
Spectrum= Container[Flux.unit*u.MeV, "dFdT"]
IntegralFlux= Container[Flux.unit*u.s*u.MeV, "dF"]

DifferentialEventRate = Container['1/(MeV*s)', "d2NdEdT"]
EventRate = Container['1/s', "dNdT"]
EventSpectrum = Container['1/MeV', "dNdE"]
EventNumber = Container['', "N"]
