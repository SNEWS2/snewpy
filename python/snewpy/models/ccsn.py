# -*- coding: utf-8 -*-
"""
The submodule ``snewpy.models.ccsn`` contains models of core-collapse supernovae
derived from the :class:`SupernovaModel` base class.

Models are initialised based on their physics parameters.
Use the ``param`` class property to view all physics parameters and their possible values:

>>> from snewpy.models.ccsn import Nakazato_2013
>>> Nakazato_2013.param
{'progenitor_mass': <Quantity [13., 20., 30., 50.] solMass>,
 'revival_time': <Quantity [  0., 100., 200., 300.] ms>,
 'metallicity': [0.02, 0.004],
 'eos': ['LS220', 'shen', 'togashi']}

For some models, not all combinations of parameters are valid. Use the ``get_param_combinations()``
class method to get a list of all valid combinations and filter it:

>>> list(params for params in Nakazato_2013.get_param_combinations() if params['eos'] != 'shen')
[{'progenitor_mass': <Quantity 30. solMass>, 'revival_time': <Quantity 0. ms>,
  'metallicity': 0.004, 'eos': 'LS220'},
 {'progenitor_mass': <Quantity 30. solMass>, 'revival_time': <Quantity 0. ms>,
  'metallicity': 0.004, 'eos': 'togashi'}]

.. _Garching Supernova Archive: https://wwwmpa.mpa-garching.mpg.de/ccsnarchive/

"""
import logging
import os
import re
import tarfile

import numpy as np
from astropy import units as u
from astropy.table import Table

from snewpy.models import ccsn_loaders as loaders

from snewpy.models.registry_model import RegistryModel, Parameter
from snewpy.models.registry_model import all_models
from textwrap import dedent


@RegistryModel()
class Fischer_2020(loaders.Fischer_2020):
    """Model based on simulations from `Fischer et al. (2020) <https://arxiv.org/abs/1804.10890>`
    """
    def __init__(self):
        self.metadata["EOS"] = "HS(DD2)"
        self.metadata["Progenitor mass"] = 18 * u.Msun
        filename='Fischer_2020.tar.gz'
        return super().__init__(filename, metadata=self.metadata)

@RegistryModel(
    progenitor_mass = [13, 20, 30, 50] * u.Msun,
    revival_time = [0, 100, 200, 300] * u.ms,
    metallicity = [0.02, 0.004],
    eos = ['LS220', 'shen', 'togashi'],

    _param_validator = lambda p: (p['revival_time'] == 0 * u.ms and p['progenitor_mass'] == 30 * u.Msun
                                  and p['metallicity'] == 0.004) or \
                                 (p['revival_time'] != 0 * u.ms and p['eos'] == 'shen'
                                  and not (p['progenitor_mass'] == 30 * u.Msun and p['metallicity'] == 0.004))
)
class Nakazato_2013(loaders.Nakazato_2013):
    """Model based on simulations from Nakazato et al., ApJ S 205:2
    (2013), ApJ 804:75 (2015), PASJ 73:639 (2021). See also http://asphwww.ph.noda.tus.ac.jp/snn/.
    """
    def __init__(self, progenitor_mass:u.Quantity, revival_time:u.Quantity, metallicity:float, eos:str):
        # Strip units for filename construction
        progenitor_mass = progenitor_mass.to(u.Msun).value
        revival_time = revival_time.to(u.ms).value

        if revival_time != 0:
            filename = f"nakazato-{eos}-z{metallicity}-t_rev{int(revival_time)}ms-s{progenitor_mass:3.1f}.fits"
        else:
            filename = f"nakazato-{eos}-BH-z{metallicity}-s{progenitor_mass:3.1f}.fits"
        #modify metadata if needed...
        #self.metadata['name']=value
        return super().__init__(filename, self.metadata)


@RegistryModel(
    progenitor_mass = [27., 9.6] * u.Msun,
    eos = ['LS220', 'SFHo']
)
class Sukhbold_2015(loaders.Sukhbold_2015):
    """Model based on simulations from Sukhbold et al., ApJ 821:38,2016. Models were shared privately by email.
    """
    def __init__(self, progenitor_mass:u.Quantity, eos:str):
        if progenitor_mass.value == 9.6:
            filename = f'sukhbold-{eos}-z{progenitor_mass.value:3.1f}.fits'
        else:
            filename = f'sukhbold-{eos}-s{progenitor_mass.value:3.1f}.fits'
        return super().__init__(filename, self.metadata)


@RegistryModel(
    progenitor_mass = [11.2, 20., 27.] * u.Msun,
    direction = [1,2,3],
    eos=['LS220'],
    _param_validator = lambda p: (p['progenitor_mass'] == 11.2 * u.Msun and p['direction'] in (1,)) or \
        (p['progenitor_mass'] == 20. * u.Msun and p['direction'] in (1,)) or \
        (p['progenitor_mass'] == 27. * u.Msun and p['direction'] in (1,2,3))
)
class Tamborra_2014(loaders.Tamborra_2014):
    """Model based on 3D simulations from `Tamborra et al., PRD 90:045032, 2014 <https://arxiv.org/abs/1406.0006>`_.
    Data files are from the `Garching Supernova Archive`_.
    """
    def __init__(self, *, progenitor_mass:u.Quantity, direction:int):
        filename = f's{progenitor_mass.value:3.1f}c_3D_dir{direction}'
        # Metadata is handled by __init__ in _GarchingArchiveModel
        return super().__init__(filename=filename, metadata=self.metadata)

@RegistryModel(
    progenitor_mass= [11.2, 27.] * u.Msun,
    eos = ['LS220']
)
class Bollig_2016(loaders.Bollig_2016):
    """Model based on simulations from `Bollig et al. (2016) <https://arxiv.org/abs/1508.00785>`_.
    Models were taken, with permission, from the `Garching Supernova Archive`_.
    """
    def __init__(self, progenitor_mass:u.Quantity, eos:str='LS220'):
        filename = f's{progenitor_mass.value:3.1f}c'
        return super().__init__(filename=filename, metadata=self.metadata)

@RegistryModel(
    progenitor_mass = [15.] * u.Msun,
    rotation = ['fast','slow','non'],
    direction = [1,2,3],
    eos=['LS220']
)
class Walk_2018(loaders.Walk_2018):
    """Model based on SASI-dominated simulations from `Walk et al.,
    PRD 98:123001, 2018 <https://arxiv.org/abs/1807.02366>`_. Data files are from
    the `Garching Supernova Archive`_.
    """
    def __init__(self, *, progenitor_mass:u.Quantity, rotation:str, direction:int):

        filename = f's{progenitor_mass.value:3.1f}c_3D_{rotation}rot_dir{direction}'
        return super().__init__(filename=filename, metadata=self.metadata)


@RegistryModel(
    progenitor_mass = [40., 75.] * u.Msun,
    direction = [1,2,3],
    eos=['LS220'],
    _param_validator = lambda p: (p['progenitor_mass'] == 75. * u.Msun and p['direction'] in (1,2)) or \
        (p['progenitor_mass'] == 40. * u.Msun and p['direction'] in (1,2,3))
)
class Walk_2019(loaders.Walk_2019):
    """Model based on SASI-dominated simulations from `Walk et al.,
    PRD 101:123013, 2019 <https://arxiv.org/abs/1910.12971>`_. Data files are
    from the `Garching Supernova Archive`_.
    """
    def __init__(self, * ,progenitor_mass:u.Quantity, direction:int):
        filename = f's{progenitor_mass.value:3.1f}c_3DBH_dir{direction}'
        return super().__init__(filename=filename, metadata=self.metadata)


@RegistryModel(
    progenitor_mass = Parameter(values=(list(range(12, 34)) +
                                 list(range(35, 61, 5)) +
                                 [70, 80, 100, 120]) * u.Msun,
                                desc_values='[12..33, 35..5..60, 70, 80, 100, 120] solMass'
                               ),
    eos = ['HShen', 'LS220']
)
class OConnor_2013(loaders.OConnor_2013):
    """Model based on the black hole formation simulation in `O'Connor & Ott (2013) <https://arxiv.org/abs/1207.1100>`_.
    """
    def __init__(self, eos:str, progenitor_mass:u.Quantity):
        # Load from Parameters
        filename = f'{eos}_timeseries.tar.gz'
        return super().__init__(filename=filename, metadata=self.metadata)

@RegistryModel(
    progenitor_mass = [40 * u.Msun],
)
class OConnor_2015(loaders.OConnor_2015):
    """Model based on the black hole formation simulation in `O'Connor (2015) <https://arxiv.org/abs/1411.7058>`_.
    """
    def __init__(self, progenitor_mass:u.Quantity):
        self.metadata["EOS"] = "LS220"
        # Filename is currently the same regardless of parameters
        filename = 'M1_neutrinos.dat'
        return super().__init__(filename, self.metadata)

@RegistryModel(
    progenitor_mass = Parameter(values=(list(range(16, 27)) + [19.89, 22.39, 30, 33]) * u.Msun,
                                desc_values = '[16..26, 19.89, 22.39, 30, 33] solMass'
                               ),
)
class Zha_2021(loaders.Zha_2021):
    """Model based on the hadron-quark phse transition models from `Zha et al. 2021 <https://arxiv.org/abs/2103.02268>`_.
    """
    def __init__(self, *, progenitor_mass:u.Quantity):
        self.metadata["EOS"] = "STOS_B145"
        filename = f's{progenitor_mass.value:g}.dat'
        return super().__init__(filename, self.metadata)

@RegistryModel(
    progenitor_mass=Parameter(np.concatenate((np.linspace(9.0, 12.75, 16),
                                              np.linspace(13, 30., 171),
                                              np.linspace(31., 33., 3),
                                              np.linspace(35, 55, 5),
                                              np.linspace(60, 80, 3),
                                              np.linspace(100, 120, 2))) << u.Msun,
                              desc_values='[9..0.25..13, 13..0.1..30, 31..33, 35..5..60, 70, 80, 100, 120] solMass'),
    turbmixing_param= Parameter([1.23, 1.25, 1.27],
                              name='turbmixing_param',
                              label='Turb. mixing param.',
                              description='Turbulent mixing parameter alpha_lambda',
                              ),
)
class Warren_2020(loaders.Warren_2020):
    """Model based on simulations from Warren et al., ApJ 898:139, 2020.
    Neutrino fluxes available at https://doi.org/10.5281/zenodo.3667908."""
    # np.arange with decimal increments can produce floating point errors
    # Though it may be more intuitive to use np.arange, these fp-errors quickly become troublesome
    def __init__(self, *, progenitor_mass, turbmixing_param):
        self.metadata["EOS"] = "SFHo"
        if progenitor_mass.value.is_integer() and progenitor_mass.value <= 30.:
            fname = f'stir_a{turbmixing_param:3.2f}/stir_multimessenger_a{turbmixing_param:3.2f}_m{progenitor_mass.value:.1f}.h5'
        else:
            fname = f'stir_a{turbmixing_param:3.2f}/stir_multimessenger_a{turbmixing_param:3.2f}_m{progenitor_mass.value:g}.h5'
        return super().__init__(fname, self.metadata)

@RegistryModel(
    rotational_velocity= [0, 1] * u.rad / u.s,
    magnetic_field_exponent= Parameter([0, 12, 13],
                                      label='B_0 Exponent',
                                      description='Exponent of magnetic field (See Eq. 46)'),
    progenitor_mass=[20*u.Msun],
    _param_validator = lambda p: (p['rotational_velocity'].value == 1 and p['magnetic_field_exponent'] in (12, 13)) or \
                               (p['rotational_velocity'].value == 0 and p['magnetic_field_exponent'] == 0)
)
class Kuroda_2020(loaders.Kuroda_2020):
    """Model based on simulations from `Kuroda et al. (2020) <https://arxiv.org/abs/2009.07733>`_."""
    def __init__(self, *, rotational_velocity, magnetic_field_exponent):
        self.metadata["EOS"] = "LS220"
        filename = f'LnuR{int(rotational_velocity.value):1d}0B{int(magnetic_field_exponent):02d}.dat'
        return super().__init__(filename, self.metadata)

@RegistryModel(
    progenitor_mass=[9, 10, 12, 13, 14, 15, 16, 19, 25, 60] * u.Msun,
)
class Fornax_2019(loaders.Fornax_2019):
    """Model based on 3D simulations from D. Vartanyan, A. Burrows, D. Radice, M.  A. Skinner and J. Dolence, MNRAS 482(1):351, 2019.
       Data available at https://www.astro.princeton.edu/~burrows/nu-emissions.3d/
    """
    def __init__(self, cache_flux=False, *, progenitor_mass):
        """
        Parameters
        ----------
        cache_flux : bool
            If true, pre-compute the flux on a fixed angular grid and store the values in a FITS file.
        """
        if progenitor_mass.value == 16:
            filename = f'lum_spec_{int(progenitor_mass.value):d}M_r250.h5'
        else:
            filename = f'lum_spec_{int(progenitor_mass.value):d}M.h5'
        return super().__init__(filename, self.metadata, cache_flux=cache_flux)

@RegistryModel(
    progenitor_mass = Parameter((list(range(12, 24)) + [25, 26, 26.99]) * u.Msun,
                                desc_values='[12..23, 25, 26, 26.99] solMass')
)
class Fornax_2021(loaders.Fornax_2021):
    """Model based on 3D simulations from D. Vartanyan, A. Burrows, D. Radice, M.  A. Skinner and J. Dolence, MNRAS 482(1):351, 2019.
       Data available at https://www.astro.princeton.edu/~burrows/nu-emissions.3d/
        """
    def __init__(self, progenitor_mass:u.Quantity):
        # Load from Parameters
        if progenitor_mass.value.is_integer():
            filename = f'lum_spec_{int(progenitor_mass.value):2d}M_r10000_dat.h5'
        else:
            filename = f'lum_spec_{progenitor_mass.value:.2f}M_r10000_dat.h5'
        return super().__init__(filename, self.metadata)


_fornax_2022_progenitors = [  '9.0',     '9.25',     '9.5',      '9.75',     '10.0',
                  '10.25',    '10.5',     '10.75',    '11.0',     '11.25',
                  '11.5',     '11.75',    '12.00.bh', '12.03.bh', '12.07.bh',
                  '12.1.bh',  '12.13',    '12.15',    '12.18.bh', '12.20.bh',
                  '12.25',    '12.33.bh', '12.40.bh', '12.45.bh', '12.50.bh',
                  '12.54.bh', '12.60.bh', '12.63',    '12.70',    '12.72.bh',
                  '12.75',    '12.80.bh', '12.85.bh', '12.90.bh', '12.93',
                  '12.97.bh', '13.00.bh', '13.05.bh', '13.11',    '13.25.bh',
                  '13.27.bh', '13.32.bh', '13.40.bh', '13.45',    '13.50.bh',
                  '13.60.bh', '13.75',    '13.82.bh', '13.90.bh', '13.96',
                  '14.01',    '14.13.bh', '14.25.bh', '14.40.bh', '14.41.bh',
                  '14.43',    '14.44.bh', '14.70.bh', '14.87.bh', '15.00.bh',
                  '15.01',    '15.04.bh', '15.05',    '15.38.bh', '16.43',
                  '16.65',    '16.99',    '17.00',    '17.07',    '17.10',
                  '17.40',    '17.48',    '17.50',    '17.51',    '17.83',
                  '18.04',    '18.05',    '18.09',    '18.10',    '18.50',
                  '19.02',    '19.56',    '19.83',    '19.99',    '20.08',
                  '20.09',    '20.18',    '20.37',    '21.00',    '21.68',
                  '22.00',    '22.30',    '22.82',    '23.00',    '23.04',
                  '23.43',    '24.00',    '25.00',    '26.00',    '26.99']

_fornax_2022_masses = [float(p.strip('.bh')) for p in _fornax_2022_progenitors] << u.Msun

@RegistryModel(progenitor_mass = _fornax_2022_masses )
class Fornax_2022(loaders.Fornax_2022):
    """Model based on 2D simulations of 100 progenitors from Tianshu Wang, David Vartanyan, Adam Burrows, and Matthew S.B. Coleman, MNRAS 517:543, 2022.
       Data available at https://www.astro.princeton.edu/~burrows/nu-emissions.2d.large/
        """
    #a mapping of mass to the progenitor
    _mass_to_progenitor = dict(zip(_fornax_2022_masses,_fornax_2022_progenitors))

    def __init__(self, progenitor_mass:u.Quantity):
        progenitor = self._mass_to_progenitor[progenitor_mass]
        self.metadata['Black hole'] = progenitor.endswith('.bh')
        filename = f'lum_spec_{progenitor}_dat.h5'
        return super().__init__(filename, self.metadata)


@RegistryModel(
               _param_validator = lambda p: \
               (p['axion_mass'] == 0 and p['axion_coupling'] == 0) or \
               (p['axion_mass'].to_value('MeV') == 100 and p['axion_coupling'].to_value('1e-10/GeV') in (2,4,10,12,14,16,20)) or \
               (p['axion_mass'].to_value('MeV') == 200 and p['axion_coupling'].to_value('1e-10/GeV') in (2,4,6,8,10,20)),

               axion_mass = Parameter(values=[0, 100, 200]<<u.MeV,
                                      description='Axion mass in units of MeV'),
               axion_coupling = Parameter(values=[0, 2, 4, 6, 8, 10, 12, 14, 16, 20]<<(1e-10/u.GeV),
                                          description='Axion-photon coupling, in units of 1e-10/GeV',
                                          precision=2 #round to 1e-12/u.GeV
                                         ),
               progenitor_mass=[20]*u.Msun
              )
class Mori_2023(loaders.Mori_2023):
    """Model based on 2D simulations with axionlike particles, K. Mori, T.  Takiwaki, K. Kotake and S. Horiuchi, Phys. Rev. D 108:063027, 2023. All models are based on the non-rotating 20 M_sun solar metallicity progenitor model from S.E. Woolsey and A. Heger, Phys. Rep. 442:269, 2007. Data from private communication.
        """
    def __init__(self, axion_mass:u.Quantity, axion_coupling:u.Quantity):
        # Make sure axion coupling is converted to units 1e-10/GeV:
        #axion_coupling = np.round(axion_coupling.to('1e-10/GeV'))

        if axion_mass == 0:
            filename = 't-prof_std.dat'
        else:
            filename = f't-prof_{axion_mass.to_value("MeV"):g}_{axion_coupling.to_value("1e-10/GeV"):g}.dat'

        # PNS mass table, from Mori+ 2023.
        mpns = { (0, 0):    1.78,
                 (100, 2):  1.77,
                 (100, 4):  1.76,
                 (100, 10): 1.77,
                 (100, 12): 1.77,
                 (100, 14): 1.77,
                 (100, 16): 1.77,
                 (100, 20): 1.74,
                 (200, 2):  1.77,
                 (200, 4):  1.76,
                 (200, 6):  1.75,
                 (200, 8):  1.74,
                 (200, 10): 1.73,
                 (200, 20): 1.62 }

        am = int(axion_mass.to_value('MeV')) if isinstance(axion_mass, u.quantity.Quantity) else int(axion_mass)
        ac = int(axion_coupling.to_value('1e-10/GeV')) if isinstance(axion_coupling, u.quantity.Quantity) else int(axion_coupling)

        pns_mass = mpns[(am,ac)]

        # Set the metadata.
        self.metadata['PNS mass'] = pns_mass*u.Msun
        return super().__init__(filename, self.metadata)

@RegistryModel(
    Bfield = ['hydro','L1','L2'],
    direction = ['average','equator','north','south'],
    grav=['A','B',None],
    rotation=[0,90,None],
    _param_validator = lambda p: (p['Bfield'] == 'hydro' and p['grav'] == None and p['rotation'] == None ) or
        (p['Bfield'] == 'L1' and p['grav'] == None and p['rotation'] in [0,90]) or
        (p['Bfield'] == 'L2' and p['rotation'] == None and p['grav'] in ['A','B'])
)
class Bugli_2021(loaders.Bugli_2021):
    """Model based on `Buggli (2021) <https://arxiv.org/abs/2105.00665>`_.
    """
    def __init__(self, *, Bfield:str, direction:str, rotation:int=None, grav:str=None):

        filename = f'{Bfield}_3d_snewpy_{direction}.dat'
        self.metadata['Progenitor mass'] = 35*u.Msun
        self.metadata['EOS'] = 'LS220'
        if Bfield=='L2':
            filename = f'{Bfield}_b12_dipdecay_3d_grav{grav}_snewpy_{direction}.dat'
        if Bfield=='L1':
            filename = f'{Bfield}_b12_3d_{rotation}deg_snewpy_{direction}.dat'
        return super().__init__(filename=filename, metadata=self.metadata)

