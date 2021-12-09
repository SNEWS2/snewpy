from __future__ import annotations
from typing import Union, Optional, List

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
            List of dimension names and their values. The order and sizes of the axes 
            should follow the dimensions of the input data
        Raises
        ------
        ValueError: Inconsistent shape
            when the dimensions of data and axes do not match
        """

        self.array = data
        self.axes= {a:u.Quantity(axes[a]) for a in axes}
        self._axnum  = {name:num for num,name in enumerate(self.axes)}
        self._axname = list(self.axes)
        self._integral = {}
        self.shape = self.array.shape
        shape_axes  = tuple(len(ax) for ax in self.axes.values())
        if self.shape!=shape_axes:
            raise ValueError(f'Inconsistent shape: array {self.shape} vs axes{shape_axes}')
            
    def __getitem__(self, args: Union[Union[int,slice], List[Union[int,slice]]]) -> Flux:
        """ Slice the flux array and produce a new flux object 
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

    def get_axis_name(self, num):
        return self._axname[num]

    def get_axis(self, axis: Union[str, int]) -> tuple[str, int]:
        """ Get the tuple of (axis name, axis number) for a given axis.
        This is a utility function.
        Parameters
        ----------
        axis: str | int
            name or number of dimension

        Returns
        -------
        tuple[str, int]
            axis name, axis number
        """
        if isinstance(axis, int):
            axname = self.get_axis_name(axis)
            axnum = axis
        else:
            axnum = self.get_axis_num(axis)
            axname = axis
        return axname, axnum
        
    def _init_integral(self, axname):
        x = self.axes[axname]
        axnum = self.get_axis_num(axname)
        cumulative = cumulative_trapezoid(self.array, x=x, axis=axnum, initial=0)
        self._integral[axname] = interp1d(x=x, y=cumulative, fill_value=0, axis = axnum, bounds_error=False)

    def squeeze(self, axis: Optional[Union[str, int]] = None) -> Flux:
        """remove a given dimension with length 1

        Parameters
        ----------
        axis: str | int | None
            name or number of dimension to be squeezed. 
            If None - remove all dimensions with len=1

        Returns
        -------
        Flux
            A flux with reduced dimension(s)
        """
        if axis is None:
            newarr  = self.array.squeeze()
            newaxes = {name:ax for name,ax in self.axes.items() if len(ax)>1}
        else:
            axname, axnum = self.get_axis(axis)
            newarr = self.array.squeeze(axnum)
            newaxes = {**self.axes}
            newaxes.pop(axname)
        return Flux(newarr, **newaxes)
    
    def sum(self, axis: Union[str, int]) -> Flux:
        """ sum flux along given axis, producing a reduced Flux object

        Parameters
        ----------
        axis: str | int
            name or number of the dimension to sum
        Returns
        -------
        Flux
            Reduced flux summed over given dimension 
            (its shape will be without that dimension)
        """
        axname, axnum = self.get_axis(axis)
        newaxes = {**self.axes}
        newaxes.pop(axname)
        newarr = self.array.sum(axnum)
        return Flux(newarr, **newaxes)

        
    def integral(self, axis: Union[str, int], limits=None) -> Flux:
        """ Integrate the flux along a specific axis `axname` within limits.

        Parameters
        ----------
        axis: str | int
            name or number of the dimension to integrate on
        limits: tuple[float, float] |  Iterable[tuple[float, float]] | None
            Integration limits, a pair of values (xMin,xMax) or list of such pairs.
            Integrate aver the whole range if limits==None

        Returns
        -------
        Flux
            The flux integrated along the given dimension with linear interpolation
        """

        axname, axnum = self.get_axis(axis)
        fint = self._integral.get(axname,self._init_integral(axname))
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
