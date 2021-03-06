# -*- coding: utf-8 -*-
"""Unit tests for the neutrino submodule.
"""

from snewpy.neutrino import Flavor
from snewpy.flavor_transformation import NoTransformation
from snewpy.models import Nakazato_2013

from astropy import units as u

import numpy as np


def test_Nakazato_vanilla():
    xform = NoTransformation()
    mfile = 'models/Nakazato_2013/nakazato-shen-z0.004-t_rev100ms-s13.0.fits'
    model = Nakazato_2013(mfile)

    assert(model.EOS == 'SHEN')
    assert(model.progenitor_mass == 13.*u.Msun)
    assert(model.revival_time == 100.*u.ms)

    assert(model.time[0] == -50*u.ms)

    tunit = model.time.unit
    Eunit = model.meanE[Flavor.NU_E].unit

    t = -50*u.ms
    assert(np.interp(t.to(tunit), model.time, model.meanE[Flavor.NU_E])*Eunit == 6.79147181061522*u.MeV)

