# -*- coding: utf-8 -*-
import unittest
import pytest

from snewpy.neutrino import MassHierarchy, MixingParameters, FourFlavorMixingParameters
import snewpy.flavor_transformation as xforms

from snewpy.flavor import ThreeFlavor, FourFlavor, TwoFlavor

from astropy import units as u
from astropy import constants as c
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
        P = xforms.NoTransformation().P(t, E)
        assert P.flavor_in == ThreeFlavor
        assert P.flavor_out == ThreeFlavor
        for f1 in ThreeFlavor:
            for f2 in ThreeFlavor:
                assert np.isclose(P[f1,f2] , 0 if f1!=f2 else 1)
                
    def test_CompleteExchange(self):
        """
        Survival probabilities for complete electron->X transformation
        """
        P = xforms.CompleteExchange().P(t, E)
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