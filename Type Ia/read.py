import tarfile
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table

scenario = 'DDT'
model = 'N100v'

def get_data(scenario, model, tbin):
    tf = tarfile.open('{}_NeutrinoFlux.tar.bz2'.format(scenario))
    dataname = '30TimeBins/NuFlux{}_TBin{:d}_NoOsc.dat'.format(model, tbin)
    datafile = tf.extractfile(dataname)

    # Process the metadata line.
    meta = datafile.readline()
    tokens = meta.split(b'@')
    t = float(tokens[1].split(b'=')[-1][:-3])

    tokens = tokens[2].split(b')(')
    dt = float(tokens[0].split(b'=')[-1][:-1])

    tokens = tokens[1].split(b'=')
    j = tokens[1].index(b'MeV')
    dE = float(tokens[1][:j])

    # Extract the data into an astropy Table.
    tab = Table.read(datafile, format='ascii.commented_header', header_start=-1)
    tab.meta['t'] = t
    tab.meta['dt'] = dt
    tab.meta['dE'] = dE

    return tab

# Plot neutrino spectra for a selection of times and for each flavor.
tbins = np.arange(1, 31, 4)
n = len(tbins)
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.plasma_r(np.linspace(0,1,n)))

fig, axes = plt.subplots(2, 3, figsize=(11,5),
                         sharex=True, sharey=True,
                         gridspec_kw={'hspace':0, 'wspace':0})
axes = axes.flatten()

for tbin in np.arange(1, 31, 4):
    tab = get_data(scenario, model, tbin)
    t = tab.meta['t']
    dt = tab.meta['dt']
    dE = tab.meta['dE']

    print(t)

    for ax, fl in zip(axes, ['NuE','NuMu','NuTau','aNuE','aNuMu','aNuTau']):
        energy = 1e3 * tab['E(GeV)']
        flux = tab[fl] / dt / dE

        ax.plot(energy, flux, label='{:.3f} s'.format(t))
        ax.set(xscale='log',
               xlim=(0.2, 50),
               yscale='log',
               ylim=(2e-3, 9e8))

axes[2].legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

axes[0].set(ylabel=r'flux [cm$^{-2}$ s$^{-1}$ MeV$^{-1}$]')
axes[3].set(ylabel=r'flux [cm$^{-2}$ s$^{-1}$ MeV$^{-1}$]')
for i in range(3,6):
    axes[i].set(xlabel=r'$E_{\nu}$ [MeV]')
    axes[i].xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%g'))

fig.tight_layout()
fig.subplots_adjust(right=0.875)

plt.show()
