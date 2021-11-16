import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from astropy import units as u

def _interpolate(x,y,z,new_x=None, new_y=None, logx=True):
    if(new_x is not None): #interpolate in T
        if(logx):
            z = interp1d(np.log(x),np.log(z)
                         ,axis=1, bounds_error=False, fill_value=0)(np.log(new_x))
            z = np.exp(z)
        else:
            z = interp1d(x,z,axis=1, bounds_error=False, fill_value=0)(new_x)
        x = new_x
    if(new_y is not None): #interpolate in E
        z = interp1d(y,z,axis=0, bounds_error=False, fill_value=0)(new_y)
        y = new_y
    return x,y,np.nan_to_num(z)

class Odrzywolek_2010(SupernovaModel):
    def __init__(self, fname):

        df = pd.read_csv(fname,delim_whitespace=True, skiprows=1, 
                     names=['step','time','Q','R','Eavg','Sigma','a','alpha','b'])
        self.a     = np.expand_dims(df.a,0)
        self.alpha = np.expand_dims(df.alpha,0)
        self.b     = np.expand_dims(df.b,0)
        self.times = df.time.values * u.s
        self.ex_ratio=0.36 #nuX/nuE ratio from Odrzywolek paper: (arXiv:astro-ph/0311012)

    def get_time(self):
        return self.times

    def get_initial_spectra(self, t, E, flavors=Flavor):
        Enu = np.expand_dims(E,1)
        self.luminosity = a*Enu**alpha*np.exp(-b*Enu) / u.MeV
        self.times = df.time.values
    
        t = t.to_value('s')
        E = E.to_value('MeV')
        ts,es,L = _interpolate(self.times,self.energies,self.luminosity,new_x=t,new_y=E)
        result  = {}
        for f in flavors:
           result[f] = L/u.
           if not f.is_electron:
            result[f]*=self.ex_ratio
        return result

class Patton_2019(SupernovaModel):
    def __init__(self, fname):
        df = pd.read_csv(fname, comment='#', delim_whitespace=True,
                         names=['time','Enu','Lnue','Lnuebar','Lnux','Lnuxbar'], 
                         usecols=range(6))
        self.df = df.set_index(['time','Enu']).unstack()
        fluxes = {}
        for flv in ['e','ebar','x','xbar']:
            table = df1[f'Lnu{flv}']
            Ts0=table.index.values*3600
            Es0=table.columns.values
            L =table.values.T
            #if flv in ['x','xbar']: L*=2
            fluxes[flv] = _interpolate(Ts0,Es0,L,new_x=Ts,new_y=Es)
        return fluxes

class Kato(SupernovaModel):
    def __init__(self, path):
    fluxes = {}
    for flv in ['nue','nueb','nux','nuxb']:
        if flv.startswith('nue'):
            pth = f'{path}/total_{flv}'
            ts,step = np.loadtxt(f'{pth}/lightcurve_{flv}_all.dat', usecols=[0,3]).T
            fname1 = f'{pth}/spe_all'
        if flv.startswith('nux'):
            pth = f'{path}/total_nux'
            ts,step = np.loadtxt(f'{pth}/lightcurve.dat', usecols=[1,0]).T
            if flv=='nux':
                fname1 = f'{pth}/spe_sum_mu_nu'
            else:
                fname1 = f'{pth}/spe_sum_mu'
        d2NdEdT = []
        for s in step:
            es,dNdE = np.loadtxt(f'{fname1}{s:05.0f}.dat').T
            d2NdEdT+=[dNdE]

        d2NdEdT = np.stack(d2NdEdT).T
        fluxes[flv] = _interpolate(ts,es,d2NdEdT,new_x=Ts, new_y=Es)
    fluxes['e'] = fluxes.pop('nue')
    fluxes['x'] = fluxes.pop('nux')
    fluxes['ebar'] = fluxes.pop('nueb')
    fluxes['xbar'] = fluxes.pop('nuxb')
    
    return fluxes

