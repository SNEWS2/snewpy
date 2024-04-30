import enum
import numpy as np
import typing

class FlavorScheme:
    def to_tex(self):
        """LaTeX-compatible string representations of flavor."""
        base = r'\nu'
        if self.is_antineutrino:
            base = r'\overline{\nu}'
        lepton = self.lepton.lower()
        if self.is_muon or self.is_tauon:
            lepton = '\\'+lepton
        return f"${base}_{{{lepton}}}$"

    @classmethod
    def _to_value(cls,key)->int:
        if isinstance(key, str):
            try:
                return cls[key].value
            except KeyError as e:
                raise KeyError(
                    f'Cannot find key "{key}" in {cls.__name__} sheme! Valid options are {list(cls)}'
                )
        elif isinstance(key, FlavorScheme):
            if not isinstance(key, cls):
                raise TypeError(f'Value {repr(key)} is not from {cls.__name__} sheme!')
            return key.value
        elif isinstance(key, typing.Iterable):
            return tuple(map(cls._to_value, key))
        return key
        
    @property
    def is_neutrino(self):
        return not self.is_antineutrino

    @property
    def is_antineutrino(self):
        return '_BAR' in self.name

    @property
    def is_electron(self):
        """Return ``True`` for ``TwoFlavor.NU_E`` and ``TwoFlavor.NU_E_BAR``."""
        return self.lepton=='E'

    @property
    def is_muon(self):
        """Return ``True`` for ``TwoFlavor.NU_E`` and ``TwoFlavor.NU_E_BAR``."""
        return self.lepton=='MU'

    @property
    def is_tauon(self):
        """Return ``True`` for ``TwoFlavor.NU_E`` and ``TwoFlavor.NU_E_BAR``."""
        return self.lepton=='TAU'

    @property
    def is_sterile(self):
        return self.lepton=='S'
        
    @property
    def lepton(self):
        return self.name.split('_')[1]
    

def _makeFlavorScheme(name:str, leptons:list[str]):
    enum_class =  enum.IntEnum(name, start=0, type=FlavorScheme,
                   names = [f'NU_{L}{BAR}' for L in leptons for BAR in ['','_BAR']]
                  )
    return enum_class

TwoFlavor = _makeFlavorScheme('TwoFlavor',['E','X'])
ThreeFlavor = _makeFlavorScheme('ThreeFlavor',['E','MU','TAU'])
FourFlavor = _makeFlavorScheme('ThreeFlavor',['E','MU','TAU','S'])


class FlavorMatrix:
    def __init__(self, 
                 array:np.ndarray,
                 flavors:FlavorScheme,
                 flavors2:FlavorScheme|None = None
                ):
                    self.array = np.asarray(array)
                    self.flavors1 = flavors
                    self.flavors2 = flavors2 or flavors
                    expected_shape = (len(self.flavors2), len(self.flavors1))
                    assert self.array.shape == expected_shape, \
                    f"Array shape {self.array.shape} mismatch expected {expected_shape}"

    def _convert_index(self, index):
        if isinstance(index, str) or (not isinstance(index,typing.Iterable)):
            index = [index]
        new_idx = [flavors._to_value(idx) for idx,flavors in zip(index, [self.flavors2,self.flavors1])]
        print(new_idx)
        return tuple(new_idx)
        
    def __getitem__(self, index):
        return self.array[self._convert_index(index)]
        
    def __setitem__(self, index, value):
        self.array[self._convert_index(index)] = value
        
    def __repr__(self):
        s = f'{self.__class__.__name__}:<{self.flavors1.__name__}->{self.flavors2.__name__}>:'
        s+='\n'+repr(self.array)
        return s
    
    def __matmul__(self, other):
        if isinstance(other, FlavorMatrix):
            data = np.tensordot(self.array, other.array, axes=[0,1])
            return FlavorMatrix(data, self.flavors1, other.flavors2)
        raise TypeError(f"Cannot multiply object of {self.__class__} by {other.__class__}")
        
    #methods for creating matrix
    @classmethod
    def zeros(cls, flavors1:FlavorScheme, flavors2:FlavorScheme|None = None):
        flavors2 = flavors2 or flavors1
        shape = (len(flavors2), len(flavors1))
        data = np.zeros(shape)
        return cls(data, flavors1, flavors2)
        
    @classmethod
    def eye(cls, flavors1:FlavorScheme, flavors2:FlavorScheme|None = None):
        flavors2 = flavors2 or flavors1
        shape = (len(flavors2), len(flavors1))
        data = np.eye(*shape)
        return cls(data, flavors1, flavors2)

    @classmethod
    def from_function(cls, flavors1:FlavorScheme, flavors2:FlavorScheme|None = None, *, function):
        flavors2 = flavors2 or flavors1
        data = [[function(f1,f2)
                 for f1 in flavors1]
                for f2 in flavors2]
        return cls(np.array(data,dtype=float), flavors1, flavors2)
    @classmethod
    def identity(cls, flavors1:FlavorScheme, flavors2:FlavorScheme|None = None):
        return cls.from_function(flavors1, flavors2, lambda f1,f2: 1.*(f1.name==f2.name))


def flavor_matrix(flavor1:FlavorScheme, flavor2:FlavorScheme=None):
    """A decorator for creating the flavor matrix from the given function"""
    flavor2 = flavor2 or flavor1
    def _decorator(func):
        return FlavorMatrix.from_function(flavor1, flavor2, function=func)
    return _decorator

#define some conversion matrices
@flavor_matrix(ThreeFlavor,TwoFlavor)
def convert_3to2(f1,f2):
    if (f1.name==f2.name):
        return 1.
    if (f1.is_neutrino==f2.is_neutrino)and(f2.lepton=='X' and f1.lepton in ['MU','TAU']):
        return 0.5
    return 0
    
@flavor_matrix(TwoFlavor,ThreeFlavor)
def convert_2to3(f1,f2):
    if (f1.name==f2.name):
        return 1.
    if (f1.is_neutrino==f2.is_neutrino)and(f1.lepton=='X' and f2.lepton in ['MU','TAU']):
        return 1
    return 0