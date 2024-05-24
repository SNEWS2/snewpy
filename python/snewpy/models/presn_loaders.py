# -*- coding: utf-8 -*-
"""
The submodule ``snewpy.models.presn_loaders`` contains classes to load pre-supernova
models from files stored on disk.
"""

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from astropy import units as u
from snewpy.models.base import SupernovaModel
from snewpy.neutrino import Flavor
from pathlib import Path

def _interp_T(t0, v0, dt=1e-3, dv=1e-10, axis=0):
    "dilog interpolation"
    lt0 = np.log(t0 + dt)
    lv0 = np.log(v0 + dv)
    lv1 = interp1d(lt0, lv0, axis=axis, fill_value=0, copy=False, bounds_error=False)
    return lambda t1: np.exp(lv1(np.log(t1 + dt)))


def _interp_E(e0, v0, axis=1):
    "linear interpolation"
    return interp1d(e0, v0, axis=axis, fill_value=0, copy=False, bounds_error=False)


def _interp_TE(times, energies, array, ax_t=1, ax_e=2):
    def _f(t, E):
        a_t = _interp_T(times, array, axis=ax_t)(t)
        a_te = _interp_E(energies, a_t, axis=ax_e)(E)
        return a_te

    return _f

class Odrzywolek_2010(SupernovaModel):
    """Set up a presupernova model, based on 
    [A. Odrzywolek and A. Heger, Acta Phys. Polon. B 41 (2010) 1611.]
    """

    def __init__(self, filename:str):
        df = pd.read_csv(
            self.request_file(filename),
            sep='\s+',
            skiprows=1,
            usecols=[1,6,7,8],
            names=["time","a","alpha","b"],
            index_col="time",
        )
        # interpolated in time
        self.df_t = _interp_T(df.index, df)
        self.factor = {}
        for f in Flavor:
            if f.is_electron:
                self.factor[f] = 1.0
            else:
                # nuX/nuE ratio from Odrzywolek paper: (arXiv:astro-ph/0311012)
                self.factor[f] = 0.19
        time = -df.index.to_numpy() << u.s
        super().__init__(time, self.metadata)

    def get_initial_spectra(self, t, E, flavors=Flavor):
        # negative t for time before SN
        t = -t.to_value("s")
        E = E.to_value("MeV")
        df = self.df_t(t)
        a, alpha, b = df.T
        Enu = np.expand_dims(E, 1)
        a = np.expand_dims(a, 0)
        alpha = np.expand_dims(alpha, 0)
        b = np.expand_dims(b, 0)
        fluence = a * Enu ** alpha * np.exp(-b * Enu) / (u.MeV * u.s)
        result = {f: fluence.T * self.factor[f] for f in flavors}
        return result


class Patton_2017(SupernovaModel):
    """Set up a presupernova model based on 
    [Kelly M. Patton et al 2017 ApJ 851 6, 
     https://doi.org/10.5281/zenodo.2598709]
    """
    def __init__(self, filename:str):
        df = pd.read_csv(
            self.request_file(filename),
            comment="#",
            sep='\s+',
            names=["time","Enu",Flavor.NU_E,Flavor.NU_E_BAR,Flavor.NU_X,Flavor.NU_X_BAR],
            usecols=range(6),
        )

        df = df.set_index(["time", "Enu"])
        times = df.index.levels[0].to_numpy()
        energies = df.index.levels[1].to_numpy()
        df = df.unstack("Enu")
        # make a 3d array
        self.array = np.stack([df[f] for f in Flavor], axis=0)
        self.interpolated = _interp_TE(
            times, energies, self.array, ax_t=1, ax_e=2
        )
        super().__init__(times << u.hour, self.metadata)

    def get_initial_spectra(self, t, E, flavors=Flavor):
        t = np.array(-t.to_value("hour"), ndmin=1)
        E = np.array(E.to_value("MeV"), ndmin=1)
        flux = self.interpolated(t, E) / (u.MeV * u.s)
        return {f: flux[f] for f in flavors}

class Kato_2017(SupernovaModel):
    """Set up a presupernova model based on 
    [Chinami Kato et al 2017 ApJ 848 48]
    """
    def __init__(self, path):
        fluxes = {}
        #reading the time steps values:
        times, step = np.loadtxt(self.request_file(f"{path}/total_nue/lightcurve_nue_all.dat"), usecols=[0, 3]).T

        file_base = {Flavor.NU_E: 'total_nue/spe_all',
                 Flavor.NU_E_BAR:'total_nueb/spe_all',
                 Flavor.NU_X:    'total_nux/spe_sum_mu_nu',
                 Flavor.NU_X_BAR:'total_nux/spe_sum_mu'
                 }
        for flv,file_base in file_base.items():
            d2NdEdT = []
            for s in step:
                energies, dNdE = np.loadtxt(
                    self.request_file(f"{path}/{file_base}{s:05.0f}.dat")
                ).T
                d2NdEdT += [dNdE]
            fluxes[flv] = np.stack(d2NdEdT)
        self.array = np.stack([fluxes[f] for f in Flavor], axis=0)
        self.interpolated = _interp_TE(
            times, energies, self.array, ax_t=1, ax_e=2
        )
        super().__init__(times << u.s, self.metadata)

    def get_initial_spectra(self, t, E, flavors=Flavor):
        t = np.array(-t.to_value("s"), ndmin=1)
        E = np.array(E.to_value("MeV"), ndmin=1)
        flux = self.interpolated(t, E) / (u.MeV * u.s)
        return {f: flux[f] for f in flavors}

class Yoshida_2016(SupernovaModel):
    """Set up a presupernova model based on 
    [Yoshida et al. (2016), PRD 93, 123012.]
    """
    def __init__(self, filename):
        with open(self.request_file(filename)) as f:
            data = []
            T = []
            while (line := f.readline()):
                if not line: break
                T += [float(line.split()[1])]
                data += [[np.loadtxt(f, max_rows=100).flatten() for i in range(4)]]
        times = np.array(T)
        energies = np.concatenate([
                np.linspace(0,10,1001)[1:],
                np.linspace(10,20,501)[1:]
            ])
        dNdEdT = np.stack(data, axis=1)
        #rearrange flavors from NU_E, NU_E_BAR,_NU_X, NU_X_BAR to current 
        indices = np.argsort([Flavor.NU_E,Flavor.NU_E_BAR, Flavor.NU_X, Flavor.NU_X_BAR]) 
        dNdEdT = dNdEdT.take(indices, axis=0)
        
        self.interpolated = _interp_TE(
            times, energies, dNdEdT, ax_t=1, ax_e=2
        )
        super().__init__(times << u.s, self.metadata)

    def get_initial_spectra(self, t, E, flavors=Flavor):
        t = np.array(-t.to_value("s"), ndmin=1)
        E = np.array(E.to_value("MeV"), ndmin=1)
        flux = self.interpolated(t, E) / (u.MeV * u.s)
        return {f: flux[f] for f in flavors}
