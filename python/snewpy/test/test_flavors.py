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
        #short notations
        assert ThreeFlavor['E'] == ThreeFlavor['e'] == ThreeFlavor.NU_E
        assert ThreeFlavor['MU'] == ThreeFlavor['mu'] == ThreeFlavor.NU_MU
        assert ThreeFlavor['MU_BAR'] == ThreeFlavor['mu_bar'] == ThreeFlavor.NU_MU_BAR
        with pytest.raises(KeyError):
            TwoFlavor['NU_MU']
        with pytest.raises(KeyError):
            ThreeFlavor['NU_X']
    @staticmethod
    def test_getitem_collective_names():
        assert ThreeFlavor['NU']==(ThreeFlavor.NU_E, ThreeFlavor.NU_MU, ThreeFlavor.NU_TAU)
        assert ThreeFlavor['NU']==ThreeFlavor['e','mu','tau']
        assert ThreeFlavor['NU_BAR']==(ThreeFlavor.NU_E_BAR, ThreeFlavor.NU_MU_BAR, ThreeFlavor.NU_TAU_BAR)
        assert ThreeFlavor['NU_BAR']==ThreeFlavor['e_bar','mu_bar','tau_bar']
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
    def test_getitem():
        m = FlavorMatrix.eye(TwoFlavor,TwoFlavor)
        assert m[TwoFlavor.NU_E, TwoFlavor.NU_E]==1
        assert m['NU_E','NU_E']==1
        assert m['NU_E','NU_X']==0
        assert np.allclose(m['NU_E'], [1,0,0,0])
        assert np.allclose(m['NU_E'], m['NU_E',:])

    @staticmethod
    def test_getitem_submatrix():
        m = FlavorMatrix.eye(TwoFlavor)
        assert np.allclose(m[['e','x'],['e','x']], [[1,0],[0,1]])
        assert np.allclose(m[:,:].array, m.array)

    @staticmethod
    def test_getitem_short():
        m = FlavorMatrix.eye(ThreeFlavor,ThreeFlavor)
        assert m['NU_E','NU_E']==m['e','e']
        assert m['NU_MU','NU_E']==m['mu','e']
        assert m['NU_E_BAR','NU_E']==m['e_bar','e']
        assert m['NU_TAU_BAR','NU_TAU']==m['tau_bar','tau']
        
    @staticmethod
    def test_setitem():
        m = FlavorMatrix.eye(TwoFlavor,TwoFlavor)
        m['NU_E']=[2,3,4,5]
        assert m['NU_E','NU_E']==2
        assert m['NU_E','NU_X']==4
        #check that nothing changed in other parts
        assert m['NU_X','NU_E']==0
        assert m['NU_X','NU_X']==1
        m['NU_E','NU_E_BAR']=123
        assert m['NU_E','NU_E_BAR']==123
        
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
            matrix = flavor>>flavor
            assert isinstance(matrix, FlavorMatrix)
            assert np.allclose(matrix.array, np.eye(len(flavor)))
    
    @staticmethod
    @pytest.mark.parametrize('flavor_in',flavor_schemes)
    @pytest.mark.parametrize('flavor_out',flavor_schemes)
    def test_conversion_matrices(flavor_in, flavor_out):
        M = flavor_in>>flavor_out
        assert M==flavor_out<<flavor_in
        assert isinstance(M, FlavorMatrix)
        assert M.flavor_in == flavor_in
        assert M.flavor_out == flavor_out

    @staticmethod
    @pytest.mark.parametrize('flavor_in',flavor_schemes)
    @pytest.mark.parametrize('flavor_out',flavor_schemes)
    def test_matrix_convert_to_flavor_method(flavor_in, flavor_out):
        M =  FlavorMatrix.eye(ThreeFlavor,ThreeFlavor)
        M1 = M.convert_to_flavor(flavor_in=flavor_in)
        assert M1.flavor_in == flavor_in
        assert M1.flavor_out == ThreeFlavor
        M1 = M.convert_to_flavor(flavor_out=flavor_out)
        assert M1.flavor_out == flavor_out
        assert M1.flavor_in == ThreeFlavor
        M1 = M.convert_to_flavor(flavor_in=flavor_in, flavor_out=flavor_out)
        assert M1.flavor_in == flavor_in
        assert M1.flavor_out == flavor_out
        #test lshift conversion methods
        assert flavor_out<<M == M.convert_to_flavor(flavor_out=flavor_out)
        assert M<<flavor_in == M.convert_to_flavor(flavor_in=flavor_in)
        assert flavor_out<<M<<flavor_in == M.convert_to_flavor(flavor_in=flavor_in, flavor_out=flavor_out)
        