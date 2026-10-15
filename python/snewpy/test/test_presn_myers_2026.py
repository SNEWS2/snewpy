"""Myers_2026: the pair files give nu_mu + nu_tau combined, so the loader's mu and tau spectra must sum to the file."""
import io
import re
import zipfile

import astropy.units as u
import numpy as np
from snewpy.flavor import ThreeFlavor
from snewpy.models import presn


def test_myers_2026_heavy_lepton_normalisation():
    model = presn.Myers_2026(progenitor_mass=15 * u.Msun, eta=0.8, bounds='core')
    with zipfile.ZipFile(model.request_file('15M_0.8_core.zip')) as zf:
        pair_names = sorted(n for n in zf.namelist() if '/pairLuminosity/' in n and n.endswith('.dat'))
        beta_names = {re.search(r'Profile(\d+)', n).group(1): n for n in zf.namelist()
                      if '/betaLuminosity/' in n and n.endswith('.dat')}
        for name in pair_names:
            text = zf.read(name).decode()
            t_cc = float(re.search(r't_cc\s*=\s*([-+0-9.eE]+)', text).group(1))
            if -100 < t_cc < -10:                              # a snapshot well inside the model's time range
                pair = np.loadtxt(io.StringIO(text), comments='#')
                beta = np.loadtxt(io.StringIO(zf.read(beta_names[re.search(r'Profile(\d+)', name).group(1)]).decode()),
                                  comments='#')
                break
    pair, beta = pair[pair[:, 0] > 0], beta[beta[:, 0] > 0]
    spectra = model.get_initial_spectra([t_cc] * u.hour, pair[:, 0] * u.MeV)
    value = lambda f: np.asarray(spectra.array[f].to_value(1 / (u.MeV * u.s))).squeeze()
    # absolute tolerance relative to the spectrum peak: the far tail (> 8 MeV, ~1e-24 of the peak) is set by the
    # interpolator's floor and log-space round-off, not by the normalisation tested here
    tol = lambda col: dict(rtol=1e-6, atol=1e-12 * col.max())

    # columns 4 and 5 of the pair file are nu_mu + nu_tau and their antineutrinos, combined
    assert np.allclose(value(ThreeFlavor.NU_MU) + value(ThreeFlavor.NU_TAU), pair[:, 3], **tol(pair[:, 3]))
    assert np.allclose(value(ThreeFlavor.NU_MU_BAR) + value(ThreeFlavor.NU_TAU_BAR), pair[:, 4], **tol(pair[:, 4]))
    assert np.allclose(value(ThreeFlavor.NU_MU), value(ThreeFlavor.NU_TAU))
    # electron flavours are beta + pair, unchanged
    assert np.allclose(value(ThreeFlavor.NU_E), beta[:, 1] + pair[:, 1], **tol(beta[:, 1] + pair[:, 1]))
    assert np.allclose(value(ThreeFlavor.NU_E_BAR), beta[:, 2] + pair[:, 2], **tol(beta[:, 2] + pair[:, 2]))
