import pytest
import numpy as np
import snewpy.flavor
from snewpy.flavor import TwoFlavor,ThreeFlavor,FourFlavor, FlavorMatrix, FlavorScheme

flavor_schemes = TwoFlavor,ThreeFlavor,FourFlavor

class TestFlavorScheme:
    @staticmethod
    def test_flavor_scheme_lengths():
        assert len(TwoFlavor)==4
        assert len(ThreeFlavor)==6
        assert len(FourFlavor)==8
        
    @staticmethod
    
    def test_getitem_string():
        assert TwoFlavor['NU_E'] == TwoFlavor.NU_E
        assert TwoFlavor['NU_X'] == TwoFlavor.NU_X
        with pytest.raises(KeyError):
            TwoFlavor['NU_MU']
            
    @staticmethod
    def test_getitem_enum():
        assert TwoFlavor[TwoFlavor.NU_E] == TwoFlavor.NU_E
        assert TwoFlavor[TwoFlavor.NU_X] == TwoFlavor.NU_X
        with pytest.raises(TypeError):
            TwoFlavor[ThreeFlavor.NU_E]
            
    @staticmethod
    def test_values_from_different_enums():
        assert TwoFlavor.NU_E==ThreeFlavor.NU_E
        assert TwoFlavor.NU_E_BAR==ThreeFlavor.NU_E_BAR
    
    @staticmethod
    def test_makeFlavorScheme():
        TestFlavor = FlavorScheme.from_lepton_names('TestFlavor',leptons=['A','B','C'])
        assert len(TestFlavor)==6
        assert [f.name for f in TestFlavor]==['NU_A','NU_A_BAR','NU_B','NU_B_BAR','NU_C','NU_C_BAR']

    @staticmethod
    def test_flavor_properties():
        f = ThreeFlavor.NU_E
        assert f.is_neutrino
        assert f.is_electron
        assert not f.is_muon
        assert not f.is_tauon
        assert f.lepton=='E'
        
        f = ThreeFlavor.NU_MU
        assert f.is_neutrino
        assert not f.is_electron
        assert f.is_muon
        assert not f.is_tauon
        assert f.lepton=='MU'

        f = ThreeFlavor.NU_E_BAR
        assert not f.is_neutrino
        assert f.is_electron
        assert not f.is_muon
        assert not f.is_tauon
        assert f.lepton=='E'

        f = ThreeFlavor.NU_MU_BAR
        assert not f.is_neutrino
        assert not f.is_electron
        assert f.is_muon
        assert not f.is_tauon
        assert f.lepton=='MU'

        f = ThreeFlavor.NU_TAU
        assert f.is_neutrino
        assert not f.is_electron
        assert not f.is_muon
        assert f.is_tauon
        assert f.lepton=='TAU'

        f = ThreeFlavor.NU_TAU_BAR
        assert not f.is_neutrino
        assert not f.is_electron
        assert not f.is_muon
        assert f.is_tauon
        assert f.lepton=='TAU'

class TestFlavorMatrix:
    @staticmethod
    def test_init_square_matrix():
        m = FlavorMatrix(array=np.ones(shape=(4,4)), flavor=TwoFlavor)
        assert m.shape == (4,4)
        assert m.flavor_in == TwoFlavor
        assert m.flavor_out == TwoFlavor
        
    @staticmethod
    def test_init_square_matrix_with_wrong_shape_raises_ValueError():
        with pytest.raises(ValueError):
            m = FlavorMatrix(array=np.ones(shape=(4,5)), flavor=TwoFlavor)
        with pytest.raises(ValueError):
            m = FlavorMatrix(array=np.ones(shape=(5,5)), flavor=TwoFlavor)
        with pytest.raises(ValueError):
            m = FlavorMatrix(array=np.ones(shape=(5,4)), flavor=TwoFlavor)

    @staticmethod
    def test_conversion_matrices_for_same_flavor_are_unity():
        for flavor in [TwoFlavor,ThreeFlavor,FourFlavor]:
            matrix = FlavorMatrix.conversion_matrix(flavor,flavor)
            print(matrix)
            assert np.allclose(matrix.array, np.eye(len(flavor)))

    @staticmethod
    @pytest.mark.parametrize('flavor_in',flavor_schemes)
    @pytest.mark.parametrize('flavor_out',flavor_schemes)
    def test_conversion_matrices(flavor_in, flavor_out):
        M = FlavorMatrix.conversion_matrix(flavor_out,flavor_in)
        assert M.flavor_in == flavor_in
        assert M.flavor_out == flavor_out