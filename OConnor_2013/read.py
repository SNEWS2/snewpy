# coding: utf-8
import tarfile
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from astropy.table import Table

eos = 'LS220'
mass = 30


def get_luminosity(mass, eos):
    # Open luminosity file.
    tf = tarfile.open('{}_timeseries.tar.gz'.format(eos))

    # Extract luminosity data.
    dataname = 's{:d}_{}_timeseries.dat'.format(mass, eos)
    datafile = tf.extractfile(dataname)
    lumdata = ascii.read(datafile, names=['t', 'Le', 'Lae', 'Lx',
                                          'Ee_avg', 'Eae_avg', 'Ex_avg',
                                          'Ee_rms', 'Eae_rms', 'Ex_rms'])
    return lumdata


def get_spectra(mass, eos):
    # Open spectra file.
    tf = tarfile.open('{}_timeseries_spectra.tar.gz'.format(eos))

    # Extract luminosity data.
    dataname = 's{:d}_{}_timeseries_spectra.dat'.format(mass, eos)
    datafile = tf.extractfile(dataname)

    t = []
    E = []
    spectra = []

    for line in datafile:
        tokens = [float(x) for x in line.strip().split()]
        if len(tokens) == 1:
            if t:
                E = _E
                spectra.append(Table([_Fe, _Fae, _Fx],
                                     names=['Fe', 'Fae', 'Fx'],
                                     meta={'t' : t[-1]}))
            t.append(tokens[0])
            _E, _Fe, _Fae, _Fx = [[] for _ in range(4)]
        elif len(tokens) == 4:
            _E.append(tokens[0])
            _Fe.append(tokens[1])
            _Fae.append(tokens[2])
            _Fx.append(tokens[3])

    spectra.append(Table([_Fe, _Fae, _Fx],
                         names=['Fe', 'Fae', 'Fx'],
                         meta={'t' : t[-1]}))

    return t, E, spectra

lumdata = get_luminosity(mass, eos)
t, E, spectra = get_spectra(mass, eos)

fig, axes = plt.subplots(1,2, figsize=(10,4), tight_layout=True)
ax = axes[0]
ax.plot(lumdata['t'], lumdata['Le'], label=r'$\mathcal{L}_{\nu_e}$')
ax.plot(lumdata['t'], lumdata['Lae'], label=r'$\mathcal{L}_{\bar{\nu}_e}$')
ax.plot(lumdata['t'], lumdata['Lx'], label=r'$\mathcal{L}_{\nu_x}$')
ax.set(xlabel='time [s]',
       ylabel='luminosity [erg s$^{-1}$]')
ax.legend()

ax = axes[1]
j = 100
print(np.shape(spectra))
spectrum = spectra[j]
ax.plot(E, spectrum['Fe'], label=r'$\nu_e$')
ax.plot(E, spectrum['Fae'], label=r'$\bar{\nu}_e$')
ax.plot(E, spectrum['Fx'], label=r'$\nu_x$')
ax.set(xlabel='energy [MeV]',
       ylabel=r'$dN/dE$ [cm$^{-2}$ s$^{-1}$ MeV$^{-1}$]')
ax.legend()

plt.show()
