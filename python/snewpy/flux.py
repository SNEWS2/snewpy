from .neutrino import Flavor
from .flavor_transformation import FlavorTransformation
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d
import numpy as np
from astropy import units as u

class Flux(object):
    """Neutrino flux container. This is a helper class to store and manipulate 
    flux tables, produced by SupernovaModels.
    """
    def __init__(self, data, **axes):
        """ 
        Parameters
        ----------
        data: np.ndarray 
            Flux values in N+1 dimensions array stored by (flavor,axis1,axis2...,axisN)
        axes: Dict[str,Iterable(float)]
            List of axes names and their values. The order and sizes of the axes 
            should follow the dimensions of the input data
        Raises
        ------
        ValueError: Inconsistent shape
            when the dimensions of data and axes do not match
        """

        self.array = data
        self.axes= {a:u.Quantity(axes[a]) for a in axes}
        self._axnum = {name:num for num,name in enumerate(self.axes)}
        self._integral = {}
        self.shape = self.array.shape
        shape_axes  = tuple(len(ax) for ax in self.axes.values())
        if self.shape!=shape_axes:
            raise ValueError(f'Inconsistent shape: array {self.shape} vs axes{shape_axes}')
            
    def __getitem__(self, args):
        """ Slice the flux array and produce a new flux object"""
        try: 
            iter(args)
        except:
            args=[args]
        args = [a if isinstance(a,slice) else slice(a,a+1) for a in args]
        array = self.array.__getitem__(tuple(args))
        newaxes = {**self.axes}
        for arg,ax in zip(args,self.axes):
            newaxes[ax] = self.axes[ax].__getitem__(arg)
        return Flux(array, **newaxes).squeeze()

    def get_axis_num(self, name):
        return self._axnum[name]
        
    def _init_integral(self, axname):
        x = self.axes[axname]
        axnum = self.get_axis_num(axname)
        cumulative = cumulative_trapezoid(self.array, x=x, axis=axnum, initial=0)
        self._integral[axname] = interp1d(x=x, y=cumulative, fill_value=0, axis = axnum, bounds_error=False)

    def squeeze(self, axname=None):
        """remove a given dimension"""
        if axname is None:
            newarr  = self.array.squeeze()
            newaxes = {name:ax for name,ax in self.axes.items() if len(ax)>1}
        else:
            axnum = self.get_axis_num(axname)
            newarr = self.array.squeeze(axnum)
            newaxes = {**self.axes}
            newaxes.pop(axname)
        return Flux(newarr, **newaxes)

    def integral(self, axname, limits=None):
        """ Integrate the flux along a specific axis `axname` within limits.
        Parameters
        ----------
        axname: str
            name of the axes to integrate on
        limits: tuple[float, float] |  Iterable[tuple[float,float]]
            Integration limits, a pair of values (xMin,xMax) or list of such pairs.

        Returns
        -------
        Flux
            The flux integrated along the given axis with linear interpolation
        """

        fint = self._integral.get(axname,self._init_integral(axname))
        axnum=self.get_axis_num(axname)
        ax = self.axes[axname]
        xmin,xmax = ax.min(), ax.max()
        if limits is None:
            limits = (xmin,xmax)
        try:
            limits = u.Quantity(limits, ndmin=2)
        except TypeError:
            limits = u.Quantity([u.Quantity(l) for l in limits], ndmin=2)

        limits = limits.clip(xmin,xmax)
        limits = limits.to(ax.unit)
        newaxes = {**self.axes}
        newaxes.pop(axname)
        newarr = fint(limits[:,1])-fint(limits[:,0])
        try:
            newarr = newarr.squeeze(axis=axnum)
            return Flux(newarr,**newaxes)
        except:
            newarr = list(newarr.swapaxes(0,axnum))
            return [Flux(arr,**newaxes) for arr in newarr]
    
    def __repr__(self):
        s = 'Flux: <'+' x '.join(f'{name}[{len(val)}]({val.min()}:{val.max()})' for name,val in self.axes.items())+'>'
        return s
