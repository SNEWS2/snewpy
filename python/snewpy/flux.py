from __future__ import annotations
from typing import Union, Optional, List

from .neutrino import Flavor
from .flavor_transformation import FlavorTransformation
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d
import numpy as np
from astropy import units as u
import itertools


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
        self.axes = {a: u.Quantity(axes[a]) for a in axes}
        self._axnum = {name: num for num, name in enumerate(self.axes)}
        self._axname = list(self.axes)
        self._integral = {}
        self.shape = self.array.shape
        shape_axes = tuple(len(ax) for ax in self.axes.values())
        if self.shape != shape_axes:
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
        except TypeError:
            args = [args]
        args = [a if isinstance(a, slice) else slice(a, a + 1) for a in args]
        array = self.array.__getitem__(tuple(args))
        newaxes = {**self.axes}
        for arg, ax in zip(args, self.axes):
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
        self._integral[axname] = interp1d(
            x=x, y=cumulative, fill_value=0, axis=axnum, bounds_error=False
        )

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
        fint = self._integral.get(axname, self._init_integral(axname))
        ax = self.axes[axname]
        xmin, xmax = ax.min(), ax.max()
        if limits is None:
            limits = (xmin, xmax)
        try:
            limits = u.Quantity(limits, ndmin=2)
        except TypeError:
            limits = u.Quantity([u.Quantity(l) for l in limits], ndmin=2)

        limits = limits.clip(xmin, xmax)
        limits = limits.to(ax.unit)
        newaxes = {**self.axes}
        newaxes.pop(axname)
        newarr = fint(limits[:, 1]) - fint(limits[:, 0])
        try:
            newarr = newarr.squeeze(axis=axnum)
            return Flux(newarr, **newaxes)
        except ValueError:
            newarr = list(newarr.swapaxes(0, axnum))
            return [Flux(arr, **newaxes) for arr in newarr]

    def __repr__(self):
        s = [
            f"{name}[{len(val)}]({val.min()}:{val.max()})"
            for name, val in self.axes.items()
        ]
        return f"Flux: <{' x '.join(s)}>"

    def _save_to_snowglobes(self, filename: str, header: str = ""):
        """save flux to GLoBES format text file.
        The flux dimensions must be <Flavor[4] x Enu[N]>
        Parameters
        ----------
        filename: str
            The path to output file
        header: str
            Optional header to be added at the beginning
        """
        energy = self.axes["Enu"]
        if isinstance(energy, u.Quantity):
            energy = energy.to_value("GeV")
        data = self.array
        if isinstance(data, u.Quantity):
            data = data.to_value("1/(MeV*cm**2*s)")
        table = {
            "E(GeV)": energy,
            "NuE": data[Flavor.NU_E],
            "NuMu": data[Flavor.NU_X],
            "NuTau": data[Flavor.NU_X],
            "aNuE": data[Flavor.NU_E_BAR],
            "aNuMu": data[Flavor.NU_X_BAR],
            "aNuTau": data[Flavor.NU_X_BAR],
        }
        header += " ".join(f"{key:>16}" for key in table)
        # Generate energy + number flux table.
        table = np.stack(list(table.values()))
        np.savetxt(filename, table.T, header=header, fmt="%17.8E", delimiter="")

    def to_snowglobes(self, filename: str, header=""):
        """Save this flux to the SNOwGLoBES format.
        Parameters
        ----------
        filename: str
            Output file name if processing flux without extra dimensions,
            or *template* output file name if extra dimensions are present.
        header: Optional[str]
            Header line(s) to be put in the beginning of the output file.
            Header can also have *template* substitutions for each bin in extra dimensions

        Returns
        -------
            list[str]
                List of produced files' names

        Flux must have dimensions ``<Flavor[4], Enu[N]>`` and
        optionally any extra dimensions ``<Extra1 [N1], Extra2 [N2]>``.
        When extra dimensions are present, this method loops over them and makes a file
        for each value in extra dimensions.
        User should provide a filename *template*, which will be filled with
        standard python ``str.format`` method, substituting variables listed below

        Template variables::

        ``{idx}``
            indices of extra dimensions bin, with ``_`` between them.
        ``{<dimension>}``
            Value of the given dimension
        ``{n_<dimension>}``
            Bin number of the given dimension

        Example::

            >>> #create flux with one extra dimension - time
            >>> flux = Flux(data=np.ones(shape=[4,10,5]),
                            Flavor=list(sorted(Flavor)),
                            Enu = np.linspace(0,100,10)*u.MeV,
                            time= np.linspace(0,10,5)*u.s)
            >>> flux
            Flux: <Flavor[4](0.0:3.0) x Enu[10](0.0 MeV:100.0 MeV) x time[5](0.0 s:10.0 s)>
            >>> flux.to_snowglobes('tmp/flux.idx{idx}.tbin{n_time}.{time}.dat')
            ['tmp/flux.idx0.tbin0.0.0_s.dat',
             'tmp/flux.idx1.tbin1.2.5_s.dat',
             'tmp/flux.idx2.tbin2.5.0_s.dat',
             'tmp/flux.idx3.tbin3.7.5_s.dat',
             'tmp/flux.idx4.tbin4.10.0_s.dat']
        """
        # prepare the iterables to loop over
        loop_axes = {ax: range(len(self.axes[ax])) for ax in self.axes}
        # do not loop over these - keep the whole slice
        loop_axes["Enu"] = [slice(None)]
        loop_axes["Flavor"] = [slice(None)]

        # get the cartesian product
        files = []
        for idx in itertools.product(*loop_axes.values()):
            # dictionaries to formatting strings
            index_dict = {"idx": "_".join([f"{i:d}" for i in idx if i != slice(None)])}
            nums_dict = {f"n_{name}": i for name, i in zip(self.axes, idx)}
            vals_dict = {
                f"{name}": str(self.axes[name][i]).replace(" ", "_")
                for name, i in zip(self.axes, idx)
            }
            for name in ["Enu", "Flavor"]:
                vals_dict.pop(name)

            fmt_dict = {**index_dict, **nums_dict, **vals_dict}

            fname = str(filename).format(**fmt_dict)
            self[idx]._save_to_snowglobes(fname, header=str(header).format(**fmt_dict))
            files += [fname]
        return files
