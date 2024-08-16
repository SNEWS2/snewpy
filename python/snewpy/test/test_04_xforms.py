# -*- coding: utf-8 -*-
import unittest
import pytest

from snewpy.neutrino import MassHierarchy, MixingParameters, ThreeFlavorMixingParameters, FourFlavorMixingParameters
import snewpy.flavor_transformation as xforms

from snewpy.flavor import ThreeFlavor, FourFlavor, TwoFlavor

from astropy import units as u
from astropy import constants as c
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time

import numpy as np
from numpy import sin,cos,exp,abs

# Dummy neutrino decay parameters; see arXiv:1910.01127.
@pytest.fixture
def params_NeutrinoDecay():
    return {
        'mass3':0.5 * u.eV/c.c**2,
        'lifetime':1 * u.day,
        'distance':10 * u.kpc
    }

@pytest.fixture
def params_QuantumDecoherence():
    return {
        'gamma3': (1e-27 * u.eV / (c.hbar * c.c)).to('1/kpc'),
        'gamma8': (1e-27 * u.eV / (c.hbar * c.c)).to('1/kpc'),
        'n': 0,
        'energy_ref': 10 * u.MeV
    }

@pytest.fixture(autouse=True)
def t():
    return np.arange(10) << u.s
    
@pytest.fixture(autouse=True)
def E():
    return np.linspace(1,100,21) << u.MeV


class TestFlavorTransformations:
    def test_NoTransformation(self):
        """
        Survival probabilities for no oscillations
        """
        P = xforms.NoTransformation().P_ff(t, E)
        assert P.flavor_in == ThreeFlavor
        assert P.flavor_out == ThreeFlavor
        for f1 in ThreeFlavor:
            for f2 in ThreeFlavor:
                assert np.isclose(P[f1,f2] , 0 if f1!=f2 else 1)
                
    def test_CompleteExchange(self):
        """
        Survival probabilities for complete electron->X transformation
        """
        P = xforms.CompleteExchange().P_ff(t, E)
        assert P.flavor_in == ThreeFlavor
        assert P.flavor_out == ThreeFlavor
        #general check
        for f1 in ThreeFlavor:
            for f2 in ThreeFlavor:
                if f1.is_neutrino!=f2.is_neutrino:
                    assert np.isclose(P[f1,f2],0)
                else:
                    assert np.isclose(P[f1,f2],0.5*(f1!=f2))

    def test_AdiabaticMSW_NMO(self):
        """
        Adiabatic MSW with normal ordering
        """
        mixpars = MixingParameters(MassHierarchy.NORMAL)
        th12, th13, th23 = mixpars.get_mixing_angles()
        
        P = xforms.AdiabaticMSW(mixpars).P(t, E)
        #convert to TwoFlavor case
        P = (TwoFlavor<<P<<TwoFlavor)

        assert np.isclose(P['E','E'] , sin(th13)**2)
        assert np.isclose(P['E','X'] , 1. - sin(th13)**2)
        assert np.isclose(P['X','X'] , 0.5*(1. + sin(th13)**2))
        assert np.isclose(P['X','E'] , 0.5*(1. - sin(th13)**2))

        assert np.isclose(P['E_bar','E_bar'] , (cos(th12)*cos(th13))**2)
        assert np.isclose(P['E_bar','X_bar'] , 1. - (cos(th12)*cos(th13))**2)
        assert np.isclose(P['X_bar','X_bar'] , 0.5*(1. + (cos(th12)*cos(th13))**2))
        assert np.isclose(P['X_bar','E_bar'] , 0.5*(1. - (cos(th12)*cos(th13))**2))

    def test_AdiabaticMSW_IMO(self):
        """
        Adiabatic MSW with inverted ordering
        """
        mixpars = MixingParameters(MassHierarchy.INVERTED)
        th12, th13, th23 = mixpars.get_mixing_angles()
        
        P = xforms.AdiabaticMSW(mixpars).P(t, E)
        #convert to TwoFlavor case
        P = (TwoFlavor<<P<<TwoFlavor)

        assert np.isclose(P['E','E'] , (sin(th12)*cos(th13))**2)
        assert np.isclose(P['E','X'] , 1. - (sin(th12)*cos(th13))**2)
        assert np.isclose(P['X','X'] , 0.5*(1. + (sin(th12)*cos(th13))**2))
        assert np.isclose(P['X','E'] , 0.5*(1. - (sin(th12)*cos(th13))**2))

        assert np.isclose(P['E_bar','E_bar'] , sin(th13)**2)
        assert np.isclose(P['E_bar','X_bar'] , 1. - sin(th13)**2)
        assert np.isclose(P['X_bar','X_bar'] , 0.5*(1. + sin(th13)**2))
        assert np.isclose(P['X_bar','E_bar'] , 0.5*(1. - sin(th13)**2))

    def test_NonAdiabaticMSWH_NMO(self):
        """
        Nonadiabatic MSW effect with normal ordering
        """
        # Test interface with default mixing angles defined in the submodule.
        mixpars = MixingParameters(MassHierarchy.NORMAL)
        th12, th13, th23 = mixpars.get_mixing_angles()
    
        P = xforms.NonAdiabaticMSWH(mixpars).P(t, E)
        #convert to TwoFlavor case
        P = (TwoFlavor<<P<<TwoFlavor)
        assert np.isclose(P['e','e'], (sin(th12)*cos(th13))**2)
        assert np.isclose(P['e','x'], 1. - (sin(th12)*cos(th13))**2)
        assert np.isclose(P['x','x'], 0.5*(1. + (sin(th12)*cos(th13))**2))
        assert np.isclose(P['x','e'], 0.5*(1. - (sin(th12)*cos(th13))**2))
    
        assert np.isclose(P['e_bar','e_bar'], (cos(th12)*cos(th13))**2)
        assert np.isclose(P['e_bar','x_bar'], 1. - (cos(th12)*cos(th13))**2)
        assert np.isclose(P['x_bar','x_bar'], 0.5*(1. + (cos(th12)*cos(th13))**2))
        assert np.isclose(P['x_bar','e_bar'], 0.5*(1. - (cos(th12)*cos(th13))**2))
    
    def test_NonAdiabaticMSWH_IMO(self):
        """
        Nonadiabatic MSW effect with inverted ordering
        """
        # Test interface with default mixing angles defined in the submodule.
        mixpars = MixingParameters(MassHierarchy.INVERTED)
        th12, th13, th23 = mixpars.get_mixing_angles()
    
        P = xforms.NonAdiabaticMSWH(mixpars).P(t, E)
        #convert to TwoFlavor case
        P = (TwoFlavor<<P<<TwoFlavor)
    
        assert np.isclose(P['e','e'], (sin(th12)*cos(th13))**2)
        assert np.isclose(P['e','x'], 1. - (sin(th12)*cos(th13))**2)
        assert np.isclose(P['x','x'], 0.5*(1. + (sin(th12)*cos(th13))**2))
        assert np.isclose(P['x','e'], 0.5*(1. - (sin(th12)*cos(th13))**2))
    
        assert np.isclose(P['e_bar','e_bar'], (cos(th12)*cos(th13))**2)
        assert np.isclose(P['e_bar','x_bar'], 1. - (cos(th12)*cos(th13))**2)
        assert np.isclose(P['x_bar','x_bar'], 0.5*(1. + (cos(th12)*cos(th13))**2))
        assert np.isclose(P['x_bar','e_bar'], 0.5*(1. - (cos(th12)*cos(th13))**2))

    def test_AdiabaticMSWes_NMO(self):
        """
        Four-neutrino adiabatic MSW with normal ordering
        """
        mixpars = FourFlavorMixingParameters(**MixingParameters('NORMAL'), theta14=1*u.deg)
        P = xforms.AdiabaticMSWes(mixpars).P(t,E)
        th12, th13, th23, th14, th24, th34 = mixpars.get_mixing_angles()
        #convert to TwoFlavor case
        P = (TwoFlavor<<P<<TwoFlavor)
        
        De1 = (cos(th12) * cos(th13) * cos(th14))**2
        De2 = (sin(th12) * cos(th13) * cos(th14))**2
        De3 = (sin(th13) * cos(th14))**2
        De4 = (sin(th14))**2
        Ds1 = (cos(th12) * cos(th13) * sin(th14))**2
        Ds2 = (sin(th12) * cos(th13) * sin(th14))**2
        Ds3 = (sin(th13) * sin(th14))**2
        Ds4 = (cos(th14))**2

        assert np.isclose(P['e','e'], De4)
        assert np.isclose(P['e','x'], De1 + De2)
        assert np.isclose(P['x','x'], (2 - De1 - De2 - Ds1 - Ds2)/2)
        assert np.isclose(P['x','e'], (1 - De4 - Ds4)/2)

        assert np.isclose(P['e_bar','e_bar'], De1)
        assert np.isclose(P['e_bar','x_bar'], De3 + De4)
        assert np.isclose(P['x_bar','x_bar'], (2 - De3 - De4 - Ds3 - Ds4)/2)
        assert np.isclose(P['x_bar','e_bar'], (1 - De1 - Ds1)/2)

    def test_AdiabaticMSWes_IMO(self):
        """
        Four-neutrino adiabatic MSW with inverted ordering
        """
        mixpars = FourFlavorMixingParameters(**MixingParameters('INVERTED'), theta14=1*u.deg)
        P = xforms.AdiabaticMSWes(mixpars).P(t,E)
        th12, th13, th23, th14, th24, th34 = mixpars.get_mixing_angles()
        #convert to TwoFlavor case
        P = (TwoFlavor<<P<<TwoFlavor)
        De1 = (cos(th12) * cos(th13) * cos(th14))**2
        De2 = (sin(th12) * cos(th13) * cos(th14))**2
        De3 = (sin(th13) * cos(th14))**2
        De4 = (sin(th14))**2
        Ds1 = (cos(th12) * cos(th13) * sin(th14))**2
        Ds2 = (sin(th12) * cos(th13) * sin(th14))**2
        Ds3 = (sin(th13) * sin(th14))**2
        Ds4 = (cos(th14))**2

        assert np.isclose(P['e','e'], De4)
        assert np.isclose(P['e','x'], De1 + De3)
        assert np.isclose(P['x','x'], (2 - De1 - De3 - Ds1 - Ds3)/2)
        assert np.isclose(P['x','e'], (1 - De4 - Ds4)/2)

        assert np.isclose(P['e_bar','e_bar'], De3)
        assert np.isclose(P['e_bar','x_bar'], De2 + De4)
        assert np.isclose(P['x_bar','x_bar'], (2 - De2 - De4 - Ds2 - Ds4)/2)
        assert np.isclose(P['x_bar','e_bar'], (1 - De3 - Ds3)/2)

    def test_NonAdiabaticMSWes_NMO(self):
        """
        Four-neutrino non-adiabatic MSW with normal ordering
        """
        mixpars = FourFlavorMixingParameters(**MixingParameters('NORMAL'), theta14=1*u.deg)
        P = xforms.NonAdiabaticMSWes(mixpars).P(t,E)
        th12, th13, th23, th14, th24, th34 = mixpars.get_mixing_angles()
        #convert to TwoFlavor case
        P = (TwoFlavor<<P<<TwoFlavor)

        De1 = (cos(th12) * cos(th13) * cos(th14))**2
        De2 = (sin(th12) * cos(th13) * cos(th14))**2
        De3 = (sin(th13) * cos(th14))**2
        De4 = (sin(th14))**2
        Ds1 = (cos(th12) * cos(th13) * sin(th14))**2
        Ds2 = (sin(th12) * cos(th13) * sin(th14))**2
        Ds3 = (sin(th13) * sin(th14))**2
        Ds4 = (cos(th14))**2

        assert np.isclose(P['e','e'], De3)
        assert np.isclose(P['e','x'], De1 + De2)
        assert np.isclose(P['x','x'], (2 - De1 - De2 - Ds1 - Ds2)/2)
        assert np.isclose(P['x','e'], (1 - De3 - Ds3)/2)

        assert np.isclose(P['e_bar','e_bar'], De1)
        assert np.isclose(P['e_bar','x_bar'], De2 + De3)
        assert np.isclose(P['x_bar','x_bar'], (2 - De2 - De3 - Ds2 - Ds3)/2)
        assert np.isclose(P['x_bar','e_bar'], (1 - De1 - Ds1)/2)

    def test_NonAdiabaticMSWes(self):
        """
        Four-neutrino non-adiabatic MSW with inverted ordering
        """
        mixpars = FourFlavorMixingParameters(**MixingParameters('INVERTED'), theta14=1*u.deg)
        P = xforms.NonAdiabaticMSWes(mixpars).P(t,E)
        th12, th13, th23, th14, th24, th34 = mixpars.get_mixing_angles()
        #convert to TwoFlavor case
        P = (TwoFlavor<<P<<TwoFlavor)

        De1 = (cos(th12) * cos(th13) * cos(th14))**2
        De2 = (sin(th12) * cos(th13) * cos(th14))**2
        De3 = (sin(th13) * cos(th14))**2
        De4 = (sin(th14))**2
        Ds1 = (cos(th12) * cos(th13) * sin(th14))**2
        Ds2 = (sin(th12) * cos(th13) * sin(th14))**2
        Ds3 = (sin(th13) * sin(th14))**2
        Ds4 = (cos(th14))**2

        assert np.isclose(P['e','e'], De2)
        assert np.isclose(P['e','x'], De1 + De3)
        assert np.isclose(P['x','x'], (2 - De1 - De3 - Ds1 - Ds3)/2)
        assert np.isclose(P['x','e'], (1 - De2 - Ds2)/2)

        assert np.isclose(P['e_bar','e_bar'], De3)
        assert np.isclose(P['e_bar','x_bar'], De1 + De2)
        assert np.isclose(P['x_bar','x_bar'], (2 - De1 - De2 - Ds1 - Ds2)/2)
        assert np.isclose(P['x_bar','e_bar'], (1 - De3 - Ds3)/2)
    
    @pytest.mark.skip(reason='someprobs in formulas?')
    def test_TwoFlavorDecoherence_NMO(self):
        """
        Two flavor decoherence with normal ordering
        """
        # Calculation that applies to 1% turbulence.
        mixpars = MixingParameters('NORMAL')
        th12, th13, th23 = mixpars.get_mixing_angles()
        
        P = xforms.TwoFlavorDecoherence(mixpars).P(t,E)
        #convert to TwoFlavor case
        P = (TwoFlavor<<P<<TwoFlavor)

        De1 = (cos(th12) * cos(th13))**2
        De2 = (sin(th12) * cos(th13))**2
        De3 = sin(th13)**2

        assert np.isclose(P['e','e'], (De2 + De3)/2)
        assert np.isclose(P['e','x'], 1 - (De2 + De3)/2)
        assert np.isclose(P['x','x'], (1 + (De2 + De3)/2)/2)
        assert np.isclose(P['x','e'], (1 - (De2 + De3)/2)/2)

        assert np.isclose(P['e_bar','e_bar'], De1)
        assert np.isclose(P['e_bar','x_bar'], 1 - De1)
        assert np.isclose(P['x_bar','x_bar'], (1 + De2)/2)
        assert np.isclose(P['x_bar','e_bar'], (1 - De2)/2)

    def test_TwoFlavorDecoherence_IMO(self):
        """
        Two flavor decoherence with inverted ordering
        """
        mixpars = MixingParameters('INVERTED')
        th12, th13, th23 = mixpars.get_mixing_angles()
        P = xforms.TwoFlavorDecoherence(mixpars).P(t,E)
        #convert to TwoFlavor case
        P = (TwoFlavor<<P<<TwoFlavor)
        
        De1 = (cos(th12) * cos(th13))**2
        De2 = (sin(th12) * cos(th13))**2
        De3 = sin(th13)**2

        assert np.isclose(P['e','e'], De2)
        assert np.isclose(P['e','x'], 1 - De2)
        assert np.isclose(P['x','x'], (1 + De2)/2)
        assert np.isclose(P['x','e'], (1 - De2)/2)

        assert np.isclose(P['e_bar','e_bar'], (De1 + De3)/2)
        assert np.isclose(P['e_bar','x_bar'], 1 - (De1 + De3)/2)
        assert np.isclose(P['x_bar','x_bar'], (1 + (De1 + De3)/2)/2)
        assert np.isclose(P['x_bar','e_bar'], (1 - (De1 + De3)/2)/2)

    def test_ThreeFlavorDecoherence(self):
        """
        Three flavor decoherence
        """
        P = xforms.ThreeFlavorDecoherence().P(t,E)
        #convert to TwoFlavor case
        P = (TwoFlavor<<P<<TwoFlavor)
        
        assert np.isclose(P['e','e'], 1./3)
        assert np.isclose(P['e','x'] , 2./3)
        assert np.isclose(P['x','x'] , 2./3)
        assert np.isclose(P['x','e'] , 1./3)
        
        assert np.isclose(P['e_bar','e_bar'] , 1./3)
        assert np.isclose(P['e_bar','x_bar'] , 2./3)
        assert np.isclose(P['x_bar','x_bar'] , 2./3)
        assert np.isclose(P['x_bar','e_bar'] , 1./3)
        
    def test_NeutrinoDecay_NMO(self, t, E):
        """
        Neutrino decay with NMO and new mixing angles 
        """
        mixpars = MixingParameters('NORMAL')
        th12, th13, th23 = mixpars.get_mixing_angles()
        mass=1*u.eV/c.c**2
        lifetime=1*u.day
        distance=10*u.kpc
        P = xforms.NeutrinoDecay(mixpars, mass=mass, tau=lifetime, dist=distance).P(t,E)
        #convert to TwoFlavor case
        P = (TwoFlavor<<P<<TwoFlavor)

        P_decay = np.exp(-mass*c.c*distance/(E*lifetime))
        De1 = (cos(th12) * cos(th13))**2
        De2 = (sin(th12) * cos(th13))**2
        De3 = sin(th13)**2
        
        # Check transition probabilities.
        prob_ee = De1*(1.-P_decay) + De3*P_decay

        assert np.allclose(P['e','e'], prob_ee)
        assert np.allclose(P['e','x'], De1 + De3)
        assert np.allclose(P['x','x'], 1 - 0.5*(De1 + De3))
        assert np.allclose(P['x','e'], 0.5*(1 - prob_ee))

        prob_exbar = De1*(1.-P_decay) + De2 + De3*P_decay

        assert np.isclose(P['e_bar','e_bar'], De3)
        assert (np.array_equal(P['e_bar','x_bar'], prob_exbar))
        assert (np.array_equal(P['x_bar','x_bar'], 1. - 0.5*prob_exbar))
        assert np.isclose(P['x_bar','e_bar'], 0.5*(1. - De3))

    def test_nudecay_nmo_default_mixing(self):
        """
        Neutrino decay with NMO and default mixing angles
        """
        # Use default mixing angles defined in the submodule.
        xform = NeutrinoDecay(mass=self.mass3, tau=self.lifetime, dist=self.distance)

        # Check transition probabilities (normal hierarchy is default).
        mixpars = MixingParameters()
        th12, th13, th23 = mixpars.get_mixing_angles()

        De1 = (cos(th12) * cos(th13))**2
        De2 = (sin(th12) * cos(th13))**2
        De3 = sin(th13)**2

        prob_ee = np.asarray([De1*(1.-exp(-xform.gamma(_E)*self.distance)) + De3*exp(-xform.gamma(_E)*self.distance) for _E in self.E])

        assert (np.array_equal(P['e','e'], prob_ee))
        assert np.isclose(P['e','x'], De1 + De3)
        assert np.isclose(P['x','x'], 1 - 0.5*(De1 + De3))
        assert (np.array_equal(P['x','e'], 0.5*(1 - prob_ee)))

        prob_exbar = np.asarray([De1*(1.-exp(-xform.gamma(_E)*self.distance)) + De2 + De3*exp(-xform.gamma(_E)*self.distance) for _E in self.E])

        assert np.isclose(P['e_bar','e_bar'], De3)
        assert (np.array_equal(P['e_bar','x_bar'], prob_exbar))
        assert (np.array_equal(P['x_bar','x_bar'], 1. - 0.5*prob_exbar))
        assert np.isclose(P['x_bar','e_bar'], 0.5*(1. - De3))

    def test_nudecay_imo(self):
        """
        Neutrino decay with IMO and new mixing angles
        """
        # Neutrino decay with IMO, overriding the default mixing angles.
        xform = NeutrinoDecay(mix_angles=(th12, th13, th23), mass=self.mass3, tau=self.lifetime, dist=self.distance, mh=MassHierarchy.INVERTED)

        De1 = (cos(th12) * cos(th13))**2
        De2 = (sin(th12) * cos(th13))**2
        De3 = sin(th13)**2

        # Check transition probabilities.
        prob_ee = np.asarray([De2*exp(-xform.gamma(_E)*self.distance) + De3*(1.-exp(-xform.gamma(_E)*self.distance)) for _E in self.E])

        assert (np.array_equal(P['e','e'], prob_ee))
        assert np.isclose(P['e','x'], De1 + De2)
        assert np.isclose(P['x','x'], 1 - 0.5*(De1 + De2))
        assert (np.array_equal(P['x','e'], 0.5*(1 - prob_ee)))

        prob_exbar = np.asarray([De1 + De2*np.exp(-xform.gamma(_E)*self.distance) +
                                 De3*(1-np.exp(-xform.gamma(_E)*self.distance)) for _E in self.E])

        assert np.isclose(P['e_bar','e_bar'], De3)
        assert (np.array_equal(P['e_bar','x_bar'], prob_exbar))
        assert (np.array_equal(P['x_bar','x_bar'], 1. - 0.5*prob_exbar))
        assert np.isclose(P['x_bar','e_bar'], 0.5*(1. - De3))

    def test_nudecay_imo_default_mixing(self):
        """
        Neutrino decay with IMO and default mixing angles
        """
        # Use default mixing angles defined in the submodule.
        xform = NeutrinoDecay(mass=self.mass3, tau=self.lifetime, dist=self.distance, mh=MassHierarchy.INVERTED)

        # Check transition probabilities.
        mixpars = MixingParameters(MassHierarchy.INVERTED)
        th12, th13, th23 = mixpars.get_mixing_angles()

        De1 = (cos(th12) * cos(th13))**2
        De2 = (sin(th12) * cos(th13))**2
        De3 = sin(th13)**2

        prob_ee = np.asarray([De2*exp(-xform.gamma(_E)*self.distance) +
                              De3*(1.-exp(-xform.gamma(_E)*self.distance)) for _E in self.E])

        assert (np.array_equal(P['e','e'], prob_ee))
        assert np.isclose(P['e','x'], De1 + De2)
        assert np.isclose(P['x','x'], 1 - 0.5*(De1 + De2))
        assert (np.array_equal(P['x','e'], 0.5*(1 - prob_ee)))

        prob_exbar = np.asarray([De1 + De2*np.exp(-xform.gamma(_E)*self.distance) +
                                 De3*(1-np.exp(-xform.gamma(_E)*self.distance)) for _E in self.E])

        assert np.isclose(P['e_bar','e_bar'], De3)
        assert (np.array_equal(P['e_bar','x_bar'], prob_exbar))
        assert (np.array_equal(P['x_bar','x_bar'], 1. - 0.5*prob_exbar))
        assert np.isclose(P['x_bar','e_bar'], 0.5*(1. - De3))

    def test_quantumdecoherence_nmo(self):
        """
        Quantum Decoherence with NMO and new mixing angles 
        """
        # Override the default mixing angles.
        xform = QuantumDecoherence(mix_angles=(th12, th13, th23), Gamma3=self.gamma3 * c.hbar * c.c, Gamma8=self.gamma8 * c.hbar * c.c, dist=self.distance, n=self.n, E0=self.energy_ref, mh=MassHierarchy.NORMAL)

        # Test computation survival and transition probabilities of mass states.
        _E = 10*u.MeV
        assert (xform.P11(_E) == 1/3 + 1/2 * np.exp(-(self.gamma3 * (_E/self.energy_ref)**self.n + self.gamma8 * (_E/self.energy_ref)**self.n / 3) * self.distance) + 1/6 * np.exp(-self.gamma8 * (_E/self.energy_ref)**self.n * self.distance))
        assert (xform.P21(_E) == 1/3 - 1/2 * np.exp(-(self.gamma3 * (_E/self.energy_ref)**self.n + self.gamma8 * (_E/self.energy_ref)**self.n / 3) * self.distance) + 1/6 * np.exp(-self.gamma8 * (_E/self.energy_ref)**self.n * self.distance))
        assert (xform.P22(_E) == 1/3 + 1/2 * np.exp(-(self.gamma3 * (_E/self.energy_ref)**self.n + self.gamma8 * (_E/self.energy_ref)**self.n / 3) * self.distance) + 1/6 * np.exp(-self.gamma8 * (_E/self.energy_ref)**self.n * self.distance))
        assert (xform.P31(_E) == 1/3 - 1/3 * np.exp(-self.gamma8 * (_E/self.energy_ref)**self.n * self.distance))
        assert (xform.P32(_E) == 1/3 - 1/3 * np.exp(-self.gamma8 * (_E/self.energy_ref)**self.n * self.distance))
        self.assertAlmostEqual(float(xform.P33(_E)), float(1/3 + 2/3 * np.exp(-self.gamma8 * (_E/self.energy_ref)**self.n * self.distance)), places=12)

        De1 = (cos(th12) * cos(th13))**2
        De2 = (sin(th12) * cos(th13))**2
        De3 = sin(th13)**2

        # Check flavor transition probabilities.
        prob_ee = np.asarray([xform.P31(_E)*De1 + xform.P32(_E)*De2 + xform.P33(_E)*De3 for _E in self.E])

        assert (np.array_equal(P['e','e'], prob_ee))
        assert (np.array_equal(P['e','x'], 1 - prob_ee))
        assert (np.array_equal(P['x','x'], 1 - 0.5*(1 - prob_ee)))
        assert (np.array_equal(P['x','e'], 0.5*(1 - prob_ee)))

        prob_eebar = np.asarray([xform.P11(_E)*De1 + xform.P21(_E)*De2 + xform.P31(_E)*De3 for _E in self.E])

        assert (np.array_equal(P['e_bar','x_bar'], 1 - prob_eebar))
        assert (np.array_equal(P['x_bar','x_bar'], 1. - 0.5*(1 - prob_eebar)))
        assert (np.array_equal(P['x_bar','e_bar'], 0.5*(1. - prob_eebar)))

    def test_quantumdecoherence_nmo_default_mixing(self):
        """
        Quantum decoherence with NMO and default mixing angles
        """
        # Use default mixing angles defined in the submodule.
        xform = QuantumDecoherence(Gamma3=self.gamma3 * c.hbar * c.c, Gamma8=self.gamma8 * c.hbar * c.c, dist=self.distance, n=self.n, E0=self.energy_ref)

        # Check transition probabilities (normal hierarchy is default).
        mixpars = MixingParameters()
        th12, th13, th23 = mixpars.get_mixing_angles()

        De1 = (cos(th12) * cos(th13))**2
        De2 = (sin(th12) * cos(th13))**2
        De3 = sin(th13)**2

        prob_ee = np.asarray([xform.P31(_E)*De1 + xform.P32(_E)*De2 + xform.P33(_E)*De3 for _E in self.E])

        assert (np.array_equal(P['e','e'], prob_ee))
        assert (np.array_equal(P['e','x'], 1 - prob_ee))
        assert (np.array_equal(P['x','x'], 1 - 0.5*(1 - prob_ee)))
        assert (np.array_equal(P['x','e'], 0.5*(1 - prob_ee)))

        prob_eebar = np.asarray([xform.P11(_E)*De1 + xform.P21(_E)*De2 + xform.P31(_E)*De3 for _E in self.E])

        assert (np.array_equal(P['e_bar','x_bar'], 1 - prob_eebar))
        assert (np.array_equal(P['x_bar','x_bar'], 1. - 0.5*(1 - prob_eebar)))
        assert (np.array_equal(P['x_bar','e_bar'], 0.5*(1. - prob_eebar)))

    def test_quantumdecoherence_imo(self):
        """
        Quantum Decoherence with IMO and new mixing angles 
        """
        # Override the default mixing angles.
        xform = QuantumDecoherence(mix_angles=(th12, th13, th23), Gamma3=self.gamma3 * c.hbar * c.c, Gamma8=self.gamma8 * c.hbar * c.c, dist=self.distance, n=self.n, E0=self.energy_ref, mh=MassHierarchy.INVERTED)

        # Test computation survival and transition probabilities of mass states.
        _E = 10*u.MeV
        assert (xform.P11(_E) == 1/3 + 1/2 * np.exp(-(self.gamma3 * (_E/self.energy_ref)**self.n + self.gamma8 * (_E/self.energy_ref)**self.n / 3) * self.distance) + 1/6 * np.exp(-self.gamma8 * (_E/self.energy_ref)**self.n * self.distance))
        assert (xform.P21(_E) == 1/3 - 1/2 * np.exp(-(self.gamma3 * (_E/self.energy_ref)**self.n + self.gamma8 * (_E/self.energy_ref)**self.n / 3) * self.distance) + 1/6 * np.exp(-self.gamma8 * (_E/self.energy_ref)**self.n * self.distance))
        assert (xform.P22(_E) == 1/3 + 1/2 * np.exp(-(self.gamma3 * (_E/self.energy_ref)**self.n + self.gamma8 * (_E/self.energy_ref)**self.n / 3) * self.distance) + 1/6 * np.exp(-self.gamma8 * (_E/self.energy_ref)**self.n * self.distance))
        assert (xform.P31(_E) == 1/3 - 1/3 * np.exp(-self.gamma8 * (_E/self.energy_ref)**self.n * self.distance))
        assert (xform.P32(_E) == 1/3 - 1/3 * np.exp(-self.gamma8 * (_E/self.energy_ref)**self.n * self.distance))
        self.assertAlmostEqual(float(xform.P33(_E)), float(1/3 + 2/3 * np.exp(-self.gamma8 * (_E/self.energy_ref)**self.n * self.distance)), places=12)

        De1 = (cos(th12) * cos(th13))**2
        De2 = (sin(th12) * cos(th13))**2
        De3 = sin(th13)**2

        # Check transition probabilities.
        prob_ee = np.asarray([xform.P22(_E)*De2 + xform.P21(_E)*De1 + xform.P32(_E)*De3 for _E in self.E])

        assert (np.array_equal(P['e','e'], prob_ee))
        assert (np.array_equal(P['e','x'], 1 - prob_ee))
        assert (np.array_equal(P['x','x'], 1 - 0.5*(1 - prob_ee)))
        assert (np.array_equal(P['x','e'], 0.5*(1 - prob_ee)))

        prob_eebar = np.asarray([xform.P31(_E)*De1 + xform.P32(_E)*De2 + xform.P33(_E)*De3 for _E in self.E])

        assert (np.array_equal(P['e_bar','x_bar'], 1 - prob_eebar))
        assert (np.array_equal(P['x_bar','x_bar'], 1. - 0.5*(1 - prob_eebar)))
        assert (np.array_equal(P['x_bar','e_bar'], 0.5*(1. - prob_eebar)))

    def test_quantumdecoherence_imo_default_mixing(self):
        """
        Quantum decoherence with IMO and default mixing angles
        """
        # Use default mixing angles defined in the submodule.
        xform = QuantumDecoherence(Gamma3=self.gamma3 * c.hbar * c.c, Gamma8=self.gamma8 * c.hbar * c.c, dist=self.distance, n=self.n, E0=self.energy_ref, mh=MassHierarchy.INVERTED)

        # Check transition probabilities.
        mixpars = MixingParameters(MassHierarchy.INVERTED)
        th12, th13, th23 = mixpars.get_mixing_angles()

        De1 = (cos(th12) * cos(th13))**2
        De2 = (sin(th12) * cos(th13))**2
        De3 = sin(th13)**2

        prob_ee = np.asarray([xform.P22(_E)*De2 + xform.P21(_E)*De1 + xform.P32(_E)*De3 for _E in self.E])

        assert (np.array_equal(P['e','e'], prob_ee))
        assert (np.array_equal(P['e','x'], 1 - prob_ee))
        assert (np.array_equal(P['x','x'], 1 - 0.5*(1 - prob_ee)))
        assert (np.array_equal(P['x','e'], 0.5*(1 - prob_ee)))

        prob_eebar = np.asarray([xform.P31(_E)*De1 + xform.P32(_E)*De2 + xform.P33(_E)*De3 for _E in self.E])

        assert (np.array_equal(P['e_bar','x_bar'], 1 - prob_eebar))
        assert (np.array_equal(P['x_bar','x_bar'], 1. - 0.5*(1 - prob_eebar)))
        assert (np.array_equal(P['x_bar','e_bar'], 0.5*(1. - prob_eebar)))

    def test_earthmatter_nmo(self):
        """Earth matter effects with normal ordering.
        """
        src = SkyCoord.from_name('Betelgeuse')
        det = EarthLocation.of_site('SuperK')
        obstime = Time('2021-5-26 14:14:00')
        altaz = src.transform_to(AltAz(obstime=obstime, location=det))

        mixing_params = ThreeFlavorMixingParameters(dm21_2=7.42e-5 * u.eV**2, dm31_2=2.517e-3 * u.eV**2,
                                                    theta12=33.45*u.deg, theta13=8.57*u.deg, theta23=49.2*u.deg,
                                                    deltaCP=197*u.deg)
        
        t = np.arange(0, 1, 0.01)<<u.s
        E = np.arange(1, 2, 1)<<u.MeV

        P = xforms.EarthMatter(mixing_params, altaz).P_fm(t, E)[:,:,0]
        m = np.asarray([
            [0.68017319, 0.29761409, 0.02221273, 0.,         0.,         0.        ],
            [0.,         0.,         0.,         0.68295904, 0.29483659, 0.02220437],
            [0.07398204, 0.36570521, 0.56031275, 0.,         0.,         0.        ],
            [0.,         0.,         0.,         0.07296648, 0.3667098,  0.56032372],
            [0.24584477, 0.3366807,  0.41747453, 0.,         0.,         0.,       ],
            [0.,         0.,         0.,         0.24407448, 0.33845361, 0.41747191]])

        assert(np.all(np.isclose(P, m)))

    def test_earthmatter_imo(self):
        """Earth matter effects with normal ordering.
        """
        src = SkyCoord.from_name('Betelgeuse')
        det = EarthLocation.of_site('SuperK')
        obstime = Time('2021-5-26 14:14:00')
        altaz = src.transform_to(AltAz(obstime=obstime, location=det))

        mixing_params = ThreeFlavorMixingParameters(dm21_2=7.42e-5 * u.eV**2, dm31_2=-2.4238e-3 * u.eV**2,
                                                    theta12=33.45*u.deg, theta13=8.6*u.deg, theta23=49.3*u.deg,
                                                    deltaCP=282*u.deg)
        
        t = np.arange(0, 1, 0.01)<<u.s
        E = np.arange(1, 2, 1)<<u.MeV

        P = xforms.EarthMatter(mixing_params, altaz).P_fm(t, E)[:,:,0]
        m = np.asarray([
             [0.68007654, 0.29757012, 0.02235334, 0.        , 0.        , 0.        ],
             [0.        , 0.        , 0.        , 0.68284785, 0.29478826, 0.02236389],
             [0.15268565, 0.28538995, 0.5619244 , 0.        , 0.        , 0.        ],
             [0.        , 0.        , 0.        , 0.15143103, 0.2866676 , 0.56190137],
             [0.16723781, 0.41703993, 0.41572226, 0.        , 0.        , 0.        ],
             [0.        , 0.        , 0.        , 0.16572112, 0.41854414, 0.41573474]])

        assert(np.all(np.isclose(P, m)))
