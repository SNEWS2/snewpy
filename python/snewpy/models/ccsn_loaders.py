# -*- coding: utf-8 -*-
"""
The submodule ``snewpy.models.ccsn_loaders`` contains classes to load core-collapse
supernova models from files stored on disk.
"""

import logging
import os
import re
import sys
import tarfile

from astropy import units as u
from astropy.table import Table, join
from astropy.io import ascii, fits
import h5py
import numpy as np
from scipy.special import gamma, lpmv

try:
    import healpy as hp
except ImportError:
    pass

from snewpy.models.base import PinchedModel, SupernovaModel
from snewpy.flavor import ThreeFlavor
from snewpy import _model_downloader

class GarchingArchiveModel(PinchedModel):
    """Subclass that reads models in the format used in the
    `Garching Supernova Archive <https://wwwmpa.mpa-garching.mpg.de/ccsnarchive/>`_."""
    def __init__(self, filename, eos='LS220', metadata={}):
        """Model Initialization.

        Parameters
        ----------
        filename : str
            Absolute or relative path to file with model data, we add nue/nuebar/nux.  This argument will be deprecated.
        eos: str
            Equation of state. Valid value is 'LS220'. This argument will be deprecated.

        Other Parameters
        ----------------
        progenitor_mass: astropy.units.Quantity
            Mass of model progenitor in units Msun. Valid values are {progenitor_mass}.
        Raises
        ------
        FileNotFoundError
            If a file for the chosen model parameters cannot be found
        ValueError
            If a combination of parameters is invalid when loading from parameters
        """
        # Read through the several ASCII files for the chosen simulation and
        # merge the data into one giant table.
        mergtab = None
        for flavor in ThreeFlavor:
            _sfx = flavor.name.replace('_', '').lower() if flavor.is_electron else "nux"
            _filename = '{}_{}_{}'.format(filename, eos, _sfx)
            _lname = 'L_{}'.format(flavor.name)
            _ename = 'E_{}'.format(flavor.name)
            _e2name = 'E2_{}'.format(flavor.name)
            _aname = 'ALPHA_{}'.format(flavor.name)

            # Open the requested filename using the model downloader.
            datafile = self.request_file(_filename)

            simtab = Table.read(datafile,
                                names=['TIME', _lname, _ename, _e2name],
                                format='ascii')
            simtab['TIME'].unit = 's'
            simtab[_lname].unit = '1e51 erg/s'
            simtab[_aname] = (2*simtab[_ename]**2 - simtab[_e2name]) / (simtab[_e2name] - simtab[_ename]**2)
            simtab[_ename].unit = 'MeV'
            del simtab[_e2name]

            if mergtab is None:
                mergtab = simtab
            else:
                mergtab = join(mergtab, simtab, keys='TIME', join_type='left')
                mergtab[_lname].fill_value = 0.
                mergtab[_ename].fill_value = 0.
                mergtab[_aname].fill_value = 0.
        simtab = mergtab.filled()
        if not metadata:
            metadata = {
                'Progenitor mass': float(os.path.basename(filename).split('s')[1].split('c')[0]) * u.Msun,
                'EOS': eos,
            }
        super().__init__(simtab, metadata)

class Nakazato_2013(PinchedModel):
    def __init__(self, filename, metadata={}):
        """Model initialization.

        Parameters
        ----------
        filename : str
            Absolute or relative path to FITS file with model data.

        Raises
        ------
        FileNotFoundError
            If a file for the chosen model parameters cannot be found
        """
        # Open the requested filename using the model downloader.
        datafile = self.request_file(filename)
        # Read FITS table using the astropy reader.
        simtab = Table.read(datafile)

        self.filename = os.path.basename(filename)
        super().__init__(simtab, metadata)


class Sukhbold_2015(Nakazato_2013):
    pass


class Tamborra_2014(GarchingArchiveModel):
    pass


class Bollig_2016(GarchingArchiveModel):
    pass


class Walk_2018(GarchingArchiveModel):
    pass


class Walk_2019(GarchingArchiveModel):
    pass


class OConnor_2013(PinchedModel):
    """Model based on the black hole formation simulation in `O'Connor & Ott (2013) <https://arxiv.org/abs/1207.1100>`_.
    """

    def __init__(self, filename, metadata={}):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to FITS file with model data.
        """
        datafile = self.request_file(filename)
        # Open luminosity file.
        with tarfile.open(datafile) as tf:
            # Extract luminosity data.
            dataname = 's{:d}_{}_timeseries.dat'.format(int(metadata['Progenitor mass'].value), metadata['EOS'])
            # Read FITS table using the astropy reader.
            simtab = ascii.read(tf.extractfile(dataname), names=['TIME', 'L_NU_E', 'L_NU_E_BAR', 'L_NU_X',
                                                                    'E_NU_E', 'E_NU_E_BAR', 'E_NU_X',
                                                                    'RMS_NU_E', 'RMS_NU_E_BAR', 'RMS_NU_X'])

        simtab['ALPHA_NU_E'] = (2.0 * simtab['E_NU_E'] ** 2 - simtab['RMS_NU_E'] ** 2) / (
                simtab['RMS_NU_E'] ** 2 - simtab['E_NU_E'] ** 2)
        simtab['ALPHA_NU_E_BAR'] = (2.0 * simtab['E_NU_E_BAR'] ** 2 - simtab['RMS_NU_E_BAR'] ** 2) / (
                simtab['RMS_NU_E_BAR'] ** 2 - simtab['E_NU_E_BAR'] ** 2)
        simtab['ALPHA_NU_X'] = (2.0 * simtab['E_NU_X'] ** 2 - simtab['RMS_NU_X'] ** 2) / (
                simtab['RMS_NU_X'] ** 2 - simtab['E_NU_X'] ** 2)

        # note, here L_NU_X is already divided by 4
        super().__init__(simtab, metadata)


class OConnor_2015(PinchedModel):
    """Model based on the black hole formation simulation in `O'Connor (2015) <https://arxiv.org/abs/1411.7058>`_.
    """

    def __init__(self, filename, metadata={}):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to FITS file with model data.
        """

        datafile = self.request_file(filename)
        simtab = Table.read(datafile,
                            names=['TIME', 'L_NU_E', 'L_NU_E_BAR', 'L_NU_X',
                                    'E_NU_E', 'E_NU_E_BAR', 'E_NU_X',
                                    'RMS_NU_E', 'RMS_NU_E_BAR', 'RMS_NU_X'],
                            format='ascii')

        header = ascii.read(simtab.meta['comments'], delimiter='=', format='no_header', names=['key', 'val'])
        tbounce = float(header['val'][0])
        simtab['TIME'] -= tbounce

        simtab['ALPHA_NU_E'] = (2.0*simtab['E_NU_E']**2 - simtab['RMS_NU_E']**2) / \
            (simtab['RMS_NU_E']**2 - simtab['E_NU_E']**2)
        simtab['ALPHA_NU_E_BAR'] = (2.0*simtab['E_NU_E_BAR']**2 - simtab['RMS_NU_E_BAR']**2) / \
            (simtab['RMS_NU_E_BAR']**2 - simtab['E_NU_E_BAR']**2)
        simtab['ALPHA_NU_X'] = (2.0*simtab['E_NU_X']**2 - simtab['RMS_NU_X']**2) / \
            (simtab['RMS_NU_X']**2 - simtab['E_NU_X']**2)

        # SYB: double-check on this factor of 4. Should be factor of 2?
        simtab['L_NU_X'] /= 4.0

        # prevent negative lums
        simtab['L_NU_E'][simtab['L_NU_E'] < 0] = 1
        simtab['L_NU_E_BAR'][simtab['L_NU_E_BAR'] < 0] = 1
        simtab['L_NU_X'][simtab['L_NU_X'] < 0] = 1

        self.filename = os.path.basename(filename)

        super().__init__(simtab, metadata)


class Zha_2021(OConnor_2015):
    pass

class Warren_2020(PinchedModel):
    def __init__(self, filename, metadata={}):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to file prefix, we add nue/nuebar/nux
        """
        # Open the requested filename using the model downloader.
        datafile = self.request_file(filename)

        # Open luminosity file.
        # Read data from HDF5 files, then store.
        f = h5py.File(datafile, 'r')

        simtab = Table()

        for i in range(len(f['nue_data']['lum'])):
            if f['sim_data']['shock_radius'][i][1] > 0.00001:
                bounce = f['sim_data']['shock_radius'][i][0]
                break

        simtab['TIME'] = f['nue_data']['lum'][:, 0] - bounce
        simtab['L_NU_E'] = f['nue_data']['lum'][:, 1] * 1e51
        simtab['L_NU_E_BAR'] = f['nuae_data']['lum'][:, 1] * 1e51
        simtab['L_NU_X'] = f['nux_data']['lum'][:, 1] * 1e51
        simtab['E_NU_E'] = f['nue_data']['avg_energy'][:, 1]
        simtab['E_NU_E_BAR'] = f['nuae_data']['avg_energy'][:, 1]
        simtab['E_NU_X'] = f['nux_data']['avg_energy'][:, 1]
        simtab['RMS_NU_E'] = f['nue_data']['rms_energy'][:, 1]
        simtab['RMS_NU_E_BAR'] = f['nuae_data']['rms_energy'][:, 1]
        simtab['RMS_NU_X'] = f['nux_data']['rms_energy'][:, 1]

        simtab['ALPHA_NU_E'] = (2.0 * simtab['E_NU_E'] ** 2 - simtab['RMS_NU_E'] ** 2) / \
            (simtab['RMS_NU_E'] ** 2 - simtab['E_NU_E'] ** 2)
        simtab['ALPHA_NU_E_BAR'] = (2.0 * simtab['E_NU_E_BAR'] ** 2 - simtab['RMS_NU_E_BAR']
                                    ** 2) / (simtab['RMS_NU_E_BAR'] ** 2 - simtab['E_NU_E_BAR'] ** 2)
        simtab['ALPHA_NU_X'] = (2.0 * simtab['E_NU_X'] ** 2 - simtab['RMS_NU_X'] ** 2) / \
            (simtab['RMS_NU_X'] ** 2 - simtab['E_NU_X'] ** 2)

        # Set model metadata.
        self.filename = os.path.basename(filename)

        super().__init__(simtab, metadata)


class Kuroda_2020(PinchedModel):
    def __init__(self, filename, metadata={}):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to file prefix, we add nue/nuebar/nux
        """

        # Open the requested filename using the model downloader.
        datafile = self.request_file(filename)
        # Read ASCII data.
        simtab = Table.read(datafile, format='ascii')

        # Get grid of model times.
        simtab['TIME'] = simtab['Tpb[ms]'] << u.ms
        for f in ["NU_E", "NU_E_BAR", "NU_X"]:
            fkey = re.sub('(E|X)_BAR', r'A\g<1>', f).lower()
            simtab[f'L_{f}'] = simtab[f'<L{fkey}>'] * 1e51 << u.erg / u.s
            simtab[f'E_{f}'] = simtab[f'<E{fkey}>'] << u.MeV
            # There is no pinch parameter so use alpha=2.0.
            simtab[f'ALPHA_{f}'] = np.full_like(simtab[f'E_{f}'].value, 2.)

        self.filename = os.path.basename(filename)

        super().__init__(simtab, metadata)


class Fornax_2019(SupernovaModel):
    def __init__(self, filename, metadata={}, cache_flux=False):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to FITS file with model data.
        cache_flux : bool
            If true, pre-compute the flux on a fixed angular grid and store the values in a FITS file.
        """
        # Set up model metadata.
        self.filename = filename
        self.metadata = metadata

        self.fluxunit = 1e50 * u.erg/(u.s*u.MeV)
        self.time = None

        # Conversion of flavor to key name in the model HDF5 file.
        self._flavorkeys = {ThreeFlavor.NU_E: 'nu0',
                            ThreeFlavor.NU_E_BAR: 'nu1',
                            ThreeFlavor.NU_MU: 'nu2',
                            ThreeFlavor.NU_MU_BAR: 'nu2',
                            ThreeFlavor.NU_TAU: 'nu2',
                            ThreeFlavor.NU_TAU_BAR: 'nu2'}

        # Read a cached flux file in FITS format or generate one.
        self.is_cached = cache_flux and 'healpy' in sys.modules
        if cache_flux and not 'healpy' in sys.modules:
            logger = logging.getLogger()
            logger.warning("No module named 'healpy'. Cannot enable caching.")

        if self.is_cached:

            self.E = {}
            self.dE = {}
            self.dLdE = {}
            self.luminosity = {}

            # Check if we're initializing on a FITS file or not.
            if filename.endswith('.fits'):
                fitsfile = filename
            else:
                fitsfile = filename.replace('h5', 'fits')

            if os.path.exists(fitsfile):
                self._read_fits(fitsfile)
                ntim, nene, npix = self.dLdE[Flavor.NU_E].shape
                self.npix = npix
                self.nside = hp.npix2nside(npix)
            else:
                with h5py.File(filename, 'r') as _h5file:
                    if self.time is None:
                        self.time = _h5file['nu0']['g0'].attrs['time'] * u.s

                    # Use a HEALPix grid with nside=4 (192 pixels) to cache the
                    # values of Y_lm(theta, phi).
                    self.nside = 4
                    self.npix = hp.nside2npix(self.nside)
                    thetac, phic = hp.pix2ang(self.nside, np.arange(self.npix))

                    Ylm = {}
                    for l in range(3):
                        Ylm[l] = {}
                        for m in range(-l, l+1):
                            Ylm[l][m] = self._real_sph_harm(l, m, thetac, phic)

                    # Store 3D tables of dL/dE for each flavor.
                    logger = logging.getLogger()
                    for flavor in ThreeFlavor:

                        key = self._flavorkeys[flavor]
                        logger.info('Caching {} for {} ({})'.format(filename, str(flavor), key))

                        self.E[flavor] = _h5file[key]['egroup'][()] * u.MeV
                        self.dE[flavor] = _h5file[key]['degroup'][()] * u.MeV

                        ntim, nene = self.E[flavor].shape
                        self.dLdE[flavor] = np.zeros((ntim, nene, self.npix), dtype=float)
                        # Loop over time bins.
                        for i in range(ntim):
                            # Loop over energy bins.
                            for j in range(nene):
                                dLdE_ij = 0.
                                # Sum over multipole moments.
                                for l in range(3):
                                    for m in range(-l, l+1):
                                        dLdE_ij += _h5file[key]['g{}'.format(j)
                                                                ]['l={} m={}'.format(l, m)][i] * Ylm[l][m]
                                self.dLdE[flavor][i][j] = dLdE_ij

                        # Integrate over energy to get L(t).
                        factor = 1. if flavor.is_electron else 0.25
                        self.dLdE[flavor] = self.dLdE[flavor] * factor * self.fluxunit
                        self.dLdE[flavor] = self.dLdE[flavor].to('erg/(s*MeV)')

                        self.luminosity[flavor] = np.sum(self.dLdE[flavor] * self.dE[flavor][:, :, np.newaxis], axis=1)

                    # Write output to FITS.
                    self._write_fits(fitsfile, overwrite=True)
        else:
            # Open the requested filename using the model downloader.
            datafile = self.request_file(filename)
            # Open HDF5 data file.
            self._h5file = h5py.File(datafile, 'r')

            # Get grid of model times in seconds.
            self.time = self._h5file['nu0']['g0'].attrs['time'] * u.s

    def _read_fits(self, filename):
        """Read cached angular data from FITS.

        Parameters
        ----------
        filename : str
            Input filename.
        """
        hdus = fits.open(filename)

        self.time = hdus['TIME'].data * u.Unit(hdus['TIME'].header['BUNIT'])

        for flavor in ThreeFlavor:
            name = str(flavor).split('.')[-1]

            ext = '{}_ENERGY'.format(name)
            self.E[flavor] = hdus[ext].data * u.Unit(hdus[ext].header['BUNIT'])

            ext = '{}_DE'.format(name)
            self.dE[flavor] = hdus[ext].data * u.Unit(hdus[ext].header['BUNIT'])

            ext = '{}_FLUX'.format(name)
            self.dLdE[flavor] = hdus[ext].data * u.Unit(hdus[ext].header['BUNIT'])
            self.dLdE[flavor] = self.dLdE[flavor].to('erg/(s*MeV)')

            self.luminosity[flavor] = np.sum(self.dLdE[flavor] * self.dE[flavor][:, :, np.newaxis], axis=1)

    def _write_fits(self, filename, overwrite=False):
        """Write angular-dependent calculated flux in FITS format.

        Parameters
        ----------
        filename : str
            Output filename.
        """
        hx = fits.HDUList()

        hdu_time = fits.PrimaryHDU(self.time.to_value('s'))
        hdu_time.header['EXTNAME'] = 'TIME'
        hdu_time.header['BUNIT'] = 'second'
        hx.append(hdu_time)

        for flavor in ThreeFlavor:
            name = str(flavor).split('.')[-1]

            hdu_E = fits.ImageHDU(self.E[flavor].to_value('MeV'))
            hdu_E.header['EXTNAME'] = '{}_ENERGY'.format(name)
            hdu_E.header['BUNIT'] = 'MeV'
            hx.append(hdu_E)

            hdu_dE = fits.ImageHDU(self.dE[flavor].to_value('MeV'))
            hdu_dE.header['EXTNAME'] = '{}_DE'.format(name)
            hdu_dE.header['BUNIT'] = 'MeV'
            hx.append(hdu_dE)

            hdu_flux = fits.ImageHDU(self.dLdE[flavor].to_value(str(self.fluxunit)))
            hdu_flux.header['EXTNAME'] = '{}_FLUX'.format(name)
            hdu_flux.header['BUNIT'] = str(self.fluxunit)
            hx.append(hdu_flux)

        hx.writeto(filename, overwrite=overwrite)

    def _fact(self, n):
        """Calculate n!.

        Parameters
        ----------
        n : int or float
            Input for computing n factorial.

        Returns
        -------
        factorial : float
            Factorial n!, computed as Gamma(n+1).
        """
        return gamma(n + 1.)

    def _real_sph_harm(self, l, m, theta, phi):
        """Compute orthonormalized real (tesseral) spherical harmonics Y_lm.

        Parameters
        ----------
        l : int
            Degree of the spherical harmonics.
        m : int
            Order of the spherical harmonics.
        theta : float or ndarray
            Input zenith angles.
        phi : float or ndarray
            Input azimuth angles.

        Returns
        -------
        Y_lm : float or ndarray
            Real-valued spherical harmonic function at theta, phi.
        """
        if m < 0:
            norm = np.sqrt((2*l + 1.)/(2*np.pi)*self._fact(l + m)/self._fact(l - m))
            return norm * lpmv(-m, l, np.cos(theta)) * np.sin(-m*phi)
        elif m == 0:
            norm = np.sqrt((2*l + 1.)/(4*np.pi))
            return norm * lpmv(0, l, np.cos(theta)) * np.ones_like(phi)
        else:
            norm = np.sqrt((2*l + 1.)/(2*np.pi)*self._fact(l - m)/self._fact(l + m))
            return norm * lpmv(m, l, np.cos(theta)) * np.cos(m*phi)

    def _get_binnedspectra(self, t, theta, phi):
        """Get binned neutrino spectrum at a particular time.

        Parameters
        ----------
        t : float or astropy.Quantity
            Time to evaluate initial and oscillated spectra.
        theta : astropy.Quantity
            Zenith angle of the spectral emission.
        phi : astropy.Quantity
            Azimuth angle of the spectral emission.

        Returns
        -------
        E : dict
            Dictionary of energy bin central values, keyed by neutrino flavor.
        dE : dict
            Dictionary of energy bin widths, keyed by neutrino flavor.
        binspec : dict
            Dictionary of binned model spectra, keyed by neutrino flavor.
        """
        E = {}
        dE = {}
        binspec = {}

        # Convert input time to a time index.
        t = t.to(self.time.unit)
        j = (np.abs(t - self.time)).argmin()

        for flavor in ThreeFlavor:
            # Cached data: read out the relevant time and angular rows.
            if self.is_cached:
                # Convert input angles to a HEALPix index.
                k = hp.ang2pix(self.nside, theta.to_value('radian'), phi.to_value('radian'))
                E[flavor] = self.E[flavor][j]
                dE[flavor] = self.dE[flavor][j]
                binspec[flavor] = self.dLdE[flavor][j, :, k]

            # Read the HDF5 input file directly and extract the spectra.
            else:
                key = self._flavorkeys[flavor]

                # Energy binning of the model for this flavor, in units of MeV.
                E[flavor] = self._h5file[key]['egroup'][j] * u.MeV
                dE[flavor] = self._h5file[key]['degroup'][j] * u.MeV

                # Storage of differential flux per energy, angle, and time.
                dLdE = np.zeros(len(E[flavor]), dtype=float)

                # Loop over energy bins.
                for ebin in range(len(E[flavor])):
                    dLdE_j = 0
                    # Sum over multipole moments.
                    for l in range(3):
                        for m in range(-l, l + 1):
                            Ylm = self._real_sph_harm(l, m, theta.to_value('radian'), phi.to_value('radian'))
                            dLdE_j += self._h5file[key]['g{}'.format(ebin)]['l={} m={}'.format(l, m)][j] * Ylm
                    dLdE[ebin] = dLdE_j

                factor = 1. if flavor.is_electron else 0.25
                binspec[flavor] = dLdE * factor * self.fluxunit
                binspec[flavor] = binspec[flavor].to('erg/(s*MeV)')

        return E, dE, binspec

    def _get_initial_spectra_dict(self, t, E, theta, phi, flavors=ThreeFlavor, interpolation='linear'):
        """Get neutrino spectra/luminosity curves before flavor transformation.

        Parameters
        ----------
        t : astropy.Quantity
            Time to evaluate initial spectra.
        E : astropy.Quantity or ndarray of astropy.Quantity
            Energies to evaluate the initial spectra.
        theta : astropy.Quantity
            Zenith angle of the spectral emission.
        phi : astropy.Quantity
            Azimuth angle of the spectral emission.
        flavors: iterable of snewpy.neutrino.Flavor
            Return spectra for these flavors only (default: all)
        interpolation : str
            Scheme to interpolate in spectra ('nearest', 'linear').

        Returns
        -------
        initialspectra : dict
            Dictionary of model spectra, keyed by neutrino flavor.
        """
        initialspectra = {}

        # Extract the binned spectra for the input t, theta, phi:
        _E, _dE, _spec = self._get_binnedspectra(t, theta, phi)

        # Avoid "division by zero" in retrieval of the spectrum.
        E[E == 0] = np.finfo(float).eps * E.unit
        logE = np.log10(E.to_value('MeV'))

        for flavor in flavors:

            # Linear interpolation in flux.
            if interpolation.lower() == 'linear':
                # Pad log(E) array with values where flux is fixed to zero.
                _logE = np.log10(_E[flavor].to_value('MeV'))
                _dlogE = np.diff(_logE)
                _logEbins = np.insert(_logE, 0, np.log10(np.finfo(float).eps * E.unit/u.MeV))
                _logEbins = np.append(_logEbins, _logE[-1] + _dlogE[-1])

                # Pad with values where flux is fixed to zero.
                _dLdE = _spec[flavor].to_value(self.fluxunit)
                _dLdE = np.insert(_dLdE, 0, 0.)
                _dLdE = np.append(_dLdE, 0.)

                initialspectra[flavor] = np.interp(logE, _logEbins, _dLdE) * self.fluxunit

            elif interpolation.lower() == 'nearest':
                _logE = np.log10(_E[flavor].to_value('MeV'))
                _dlogE = np.diff(_logE)[0]
                _logEbins = _logE - _dlogE
                _logEbins = np.concatenate((_logEbins, [_logE[-1] + _dlogE]))
                _Ebins = 10**_logEbins

                idx = np.searchsorted(_Ebins, E) - 1
                select = (idx > 0) & (idx < len(_E[flavor]))

                _dLdE = np.zeros(len(E))
                _dLdE[np.where(select)] = np.asarray([_spec[flavor][i].to_value(self.fluxunit) for i in idx[select]])
                initialspectra[flavor] = _dLdE * self.fluxunit

            else:
                raise ValueError('Unrecognized interpolation type "{}"'.format(interpolation))

        return initialspectra


class Fornax_2021(SupernovaModel):
    def __init__(self, filename, metadata={}):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to HDF5 file with model data.
        """

        # Open the requested filename using the model downloader.
        datafile = self.request_file(filename)
        # Set up model metadata.
        self.progenitor_mass = float(filename.split('/')[-1].split('_')[2][:-1]) * u.Msun
        self.metadata = metadata
        # Open HDF5 data file.
        _h5file = h5py.File(datafile, 'r')

        self.time = _h5file['nu0'].attrs['time'] * u.s

        self.luminosity = {}
        self._E = {}
        self._dLdE = {}
        for flavor in ThreeFlavor:
            # Convert flavor to key name in the model HDF5 file
            key = {ThreeFlavor.NU_E: 'nu0',
                   ThreeFlavor.NU_E_BAR: 'nu1',
                   ThreeFlavor.NU_MU: 'nu2',
                   ThreeFlavor.NU_MU_BAR: 'nu2',
                   ThreeFlavor.NU_TAU: 'nu2',
                   ThreeFlavor.NU_TAU_BAR: 'nu2'}[flavor]

            self._E[flavor] = np.asarray(_h5file[key]['egroup'])
            self._dLdE[flavor] = {f"g{i}": np.asarray(_h5file[key][f'g{i}']) for i in range(12)}

            # Compute luminosity by integrating over model energy bins.
            dE = np.asarray(_h5file[key]['degroup'])
            n = len(dE[0])
            dLdE = np.zeros((len(self.time), n), dtype=float)
            for i in range(n):
                dLdE[:, i] = self._dLdE[flavor][f"g{i}"]

            # Note factor of 0.25 in nu_x and nu_x_bar.
            factor = 1. if flavor.is_electron else 0.25
            self.luminosity[flavor] = np.sum(dLdE*dE, axis=1) * factor * 1e50 * u.erg/u.s

    def _get_initial_spectra_dict(self, t, E, flavors=ThreeFlavor, interpolation='linear'):
        """Get neutrino spectra/luminosity curves after oscillation.

        Parameters
        ----------
        t : astropy.Quantity
            Time to evaluate initial spectra.
        E : astropy.Quantity or ndarray of astropy.Quantity
            Energies to evaluate the initial spectra.
        flavors: iterable of snewpy.neutrino.Flavor
            Return spectra for these flavors only (default: all)
        interpolation : str
            Scheme to interpolate in spectra ('nearest', 'linear').

        Returns
        -------
        initialspectra : dict
            Dictionary of model spectra, keyed by neutrino flavor.
        """
        initialspectra = {}

        # Avoid "division by zero" in retrieval of the spectrum.
        E[E == 0] = np.finfo(float).eps * E.unit
        logE = np.log10(E.to_value('MeV'))

        # Make sure the input time uses the same units as the model time grid.
        # Convert input time to a time index.
        t = u.Quantity(t.to(self.time.unit), ndmin=1)
        j = np.array(list(np.abs(_t - self.time).argmin() for _t in t))

        for flavor in flavors:
            # Energy bin centers (in MeV)
            _E = self._E[flavor][j]
            _logE = np.log10(_E)
            _dlogE = np.diff(_logE)

            # Model flavors (internally) are nu_e, nu_e_bar, and nu_x, which stands
            # for nu_mu(_bar) and nu_tau(_bar), making the flux 4x higher than nu_e and nu_e_bar.
            factor = 1. if flavor.is_electron else 0.25

            # Linear interpolation in flux.
            if interpolation.lower() == 'linear':
                # Pad log(E) array with values where flux is fixed to zero.
                _logEbins = np.insert(_logE, 0, np.log10(np.finfo(float).eps * E.unit/u.MeV), axis=1)
                _logEbins = np.append(_logEbins, np.expand_dims(_logE[:,-1] + _dlogE[:,-1], 1), axis=1)

                # Luminosity spectrum _dLdE is in units of 1e50 erg/s/MeV.
                # Pad with values where flux is fixed to zero, then divide by E to get number luminosity
                _dNLdE = np.asarray([np.zeros(j.shape)] + [self._dLdE[flavor]['g{}'.format(i)][j] for i in range(12)] + [np.zeros(j.shape)]).T
                interp_values = np.array([np.interp(logE, __logEbins, __dNLdE)
                                          for __logEbins, __dNLdE in zip(_logEbins, _dNLdE)])
                initialspectra[flavor] = (interp_values / E * factor * 1e50 * u.erg/u.s/u.MeV).to('1 / (erg s)')

            elif interpolation.lower() == 'nearest':
                # Find edges of energy bins and identify which energy bin (each entry of) E falls into
                _logEbinEdges = _logE - _dlogE[0,0] / 2
                _logEbinEdges = np.append(_logEbinEdges, np.expand_dims(_logE[:,-1] + _dlogE[:,-1]/2, 1), axis=1)
                _EbinEdges = 10**_logEbinEdges
                idx = np.array([np.searchsorted(edges, E) - 1 for edges in _EbinEdges])
                select = np.array([(_idx > 0) & (_idx < len(__E)) for _idx, __E in zip(idx, _E)])

                # Divide luminosity spectrum by energy at bin center to get number luminosity spectrum
                _dNLdE = np.zeros([len(j), len(np.atleast_1d(E))])
                for i in range(len(j)):
                    _dNLdE[i][np.where(select[i])] = np.asarray([self._dLdE[flavor]['g{}'.format(ebin_idx)][j[i]] / _E[i][ebin_idx]
                                                                 for ebin_idx in idx[i][select[i]]])
                initialspectra[flavor] = ((_dNLdE << 1/u.MeV) * factor * 1e50 * u.erg/u.s/u.MeV).to('1 / (erg s)')

            else:
                raise ValueError('Unrecognized interpolation type "{}"'.format(interpolation))

        return initialspectra


class Fornax_2022(Fornax_2021):
    def __init__(self, filename, metadata={}):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to HDF5 file with model data.
        """
        # Open the requested filename using the model downloader.
        datafile = self.request_file(filename)
        # Set up model metadata.
        self.progenitor = os.path.splitext(os.path.basename(filename))[0].split('_')[2]
        self.progenitor_mass = float(self.progenitor[:-3])*u.Msun if self.progenitor.endswith('bh') else float(self.progenitor)*u.Msun

        self.metadata = metadata

        # Open HDF5 data file.
        _h5file = h5py.File(datafile, 'r')

        self.metadata['PNS mass'] = _h5file.attrs['Mpns'] * u.Msun
        self.time = _h5file['nu0'].attrs['time'] * u.s

        self.luminosity = {}
        self._E = {}
        self._dLdE = {}
        for flavor in ThreeFlavor:
            # Convert flavor to key name in the model HDF5 file
            key = {ThreeFlavor.NU_E: 'nu0',
                   ThreeFlavor.NU_E_BAR: 'nu1',
                   ThreeFlavor.NU_MU: 'nu2',
                   ThreeFlavor.NU_MU_BAR: 'nu2',
                   ThreeFlavor.NU_TAU: 'nu2',
                   ThreeFlavor.NU_TAU_BAR: 'nu2'}[flavor]

            self._E[flavor] = np.asarray(_h5file[key]['egroup'])
            self._dLdE[flavor] = {f"g{i}": np.asarray(_h5file[key][f'g{i}']) for i in range(12)}

            # Compute luminosity by integrating over model energy bins.
            dE = np.asarray(_h5file[key]['degroup'])
            n = len(dE[0])
            dLdE = np.zeros((len(self.time), n), dtype=float)
            for i in range(n):
                dLdE[:, i] = self._dLdE[flavor][f"g{i}"]

            # Note factor of 0.25 in nu_x and nu_x_bar.
            factor = 1. if flavor.is_electron else 0.25
            self.luminosity[flavor] = np.sum(dLdE*dE, axis=1) * factor * 1e50 * u.erg/u.s


class Mori_2023(PinchedModel):
    def __init__(self, filename, metadata={}):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to file prefix.
        """
        # Open the requested filename using the model downloader.
        datafile = self.request_file(filename)

        self.metadata = metadata

        # Read ASCII data.
        simtab = Table.read(datafile, format='ascii')

        # Remove the first table row, which appears to have zero input.
        simtab = simtab[simtab['1:t_sim[s]'] > 0]

        # Get grid of model times.
        simtab['TIME'] = simtab['2:t_pb[s]'] << u.s
        for j, (f, fkey) in enumerate(zip(["NU_E", "NU_E_BAR", "NU_X"], 'ebx')):
            simtab[f'L_{f}'] = simtab[f'{6+j}:Le{fkey}[e/s]'] << u.erg / u.s
            # Compute the pinch parameter from E_rms and E_avg
            # <E^2> / <E>^2 = (2+a)/(1+a), where
            # E_rms^2 = <E^2> - <E>^2.
            Eavg = simtab[f'{9+j}:Em{fkey}[MeV]']
            Erms = simtab[f'{12+j}:Er{fkey}[MeV]']
            E2 = Erms**2 + Eavg**2
            x = E2 / Eavg**2
            alpha = (2-x) / (x-1)

            simtab[f'E_{f}'] = Eavg << u.MeV
            simtab[f'E2_{f}'] = E2 << u.MeV**2
            simtab[f'ALPHA_{f}'] = alpha

#            simtab[f'E_{f.name}'] = simtab[f'{9+j}:Em{fkey}[MeV]'] << u.MeV
#            Erms = simtab[f'{12+j}:Er{fkey}[MeV]'] * u.MeV
#
#            # Compute the pinch parameter from E_rms and E_avg
#            simtab[f'E2_{f.name}'] = Erms**2 + simtab[f'E_{f.name}']**2
#            x = simtab[f'E2_{f.name}'] / simtab[f'E_{f.name}']**2
#            simtab[f'ALPHA_{f.name}'] = (2-x) / (x-1)

        self.filename = os.path.basename(filename)

        super().__init__(simtab, metadata)


class Bugli_2021(PinchedModel):
    """Model based on `Buggli (2021) <https://arxiv.org/abs/2105.00665>`_.
    """

    def __init__(self, filename, metadata={}):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to FITS file with model data.
        """

        datafile = self.request_file(filename)
        simtab = Table.read(datafile,
                            names=['TIME', 'L_NU_E', 'L_NU_E_BAR', 'L_NU_X',
                                    'E_NU_E', 'E_NU_E_BAR', 'E_NU_X',
                                    'RMS_NU_E', 'RMS_NU_E_BAR', 'RMS_NU_X'],
                            format='ascii')

        simtab['ALPHA_NU_E'] = (2.0*simtab['E_NU_E']**2 - simtab['RMS_NU_E']**2) / \
            (simtab['RMS_NU_E']**2 - simtab['E_NU_E']**2)
        simtab['ALPHA_NU_E_BAR'] = (2.0*simtab['E_NU_E_BAR']**2 - simtab['RMS_NU_E_BAR']**2) / \
            (simtab['RMS_NU_E_BAR']**2 - simtab['E_NU_E_BAR']**2)
        simtab['ALPHA_NU_X'] = (2.0*simtab['E_NU_X']**2 - simtab['RMS_NU_X']**2) / \
            (simtab['RMS_NU_X']**2 - simtab['E_NU_X']**2)

        self.filename = os.path.basename(filename)

        super().__init__(simtab, metadata)


class Fischer_2020(PinchedModel):
    def __init__(self, filename, metadata={}):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to file
        """
        # Open the requested filename using the model downloader.
        datafile = self.request_file(filename)
        self.metadata = metadata

        # Open the requested filename using the model downloader.
        # datafile = _model_downloader.get_model_data(self.__class__.__name__, filename)
        # self.filename = os.path.basename(filename)

        simtab = Table()

        tf = tarfile.open(datafile)

        # Open luminosity file
        with tf.extractfile("luminosity.dat") as Lfile:
             Ldata = np.genfromtxt(Lfile, skip_header=2)

             simtab['TIME'] = Ldata[:, 0]

             simtab['L_NU_E'] = Ldata[:, 1]
             simtab['L_NU_E_BAR'] = Ldata[:, 2]
             simtab['L_NU_X'] = Ldata[:, 3]
             simtab['L_NU_X_BAR'] = Ldata[:, 4]

             Lfile.close()

        # Open mean energy file
        with tf.extractfile("menergy.dat") as Efile:
            Edata = np.genfromtxt(Efile, skip_header=2)

            simtab['E_NU_E'] = Edata[:, 1] << u.MeV
            simtab['E_NU_E_BAR'] = Edata[:, 2] << u.MeV
            simtab['E_NU_X'] = Edata[:, 3] << u.MeV
            simtab['E_NU_X_BAR'] = Edata[:, 4] << u.MeV

            Efile.close()

        # Open rms energy file
        with tf.extractfile("rmsenergy.dat") as RMSEfile:
            RMSEdata = np.genfromtxt(RMSEfile, skip_header=2)

            simtab['RMS_NU_E'] = RMSEdata[:, 1] << u.MeV
            simtab['RMS_NU_E_BAR'] = RMSEdata[:, 2] << u.MeV
            simtab['RMS_NU_X'] = RMSEdata[:, 3] << u.MeV
            simtab['RMS_NU_X_BAR'] = RMSEdata[:, 4] << u.MeV

            RMSEfile.close()

        simtab['ALPHA_NU_E'] = (2.0 * simtab['E_NU_E'] ** 2 - simtab['RMS_NU_E'] ** 2) / \
            (simtab['RMS_NU_E'] ** 2 - simtab['E_NU_E'] ** 2)
        simtab['ALPHA_NU_E_BAR'] = (2.0 * simtab['E_NU_E_BAR'] ** 2 - simtab['RMS_NU_E_BAR'] ** 2) / \
            (simtab['RMS_NU_E_BAR'] ** 2 - simtab['E_NU_E_BAR'] ** 2)
        simtab['ALPHA_NU_X'] = (2.0 * simtab['E_NU_X'] ** 2 - simtab['RMS_NU_X'] ** 2) / \
            (simtab['RMS_NU_X'] ** 2 - simtab['E_NU_X'] ** 2)
        simtab['ALPHA_NU_X_BAR'] = (2.0 * simtab['E_NU_X_BAR'] ** 2 - simtab['RMS_NU_X_BAR'] ** 2) / \
            (simtab['RMS_NU_X_BAR'] ** 2 - simtab['E_NU_X_BAR'] ** 2)

        tf.close()

        super().__init__(simtab, metadata)
