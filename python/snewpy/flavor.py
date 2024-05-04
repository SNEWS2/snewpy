import enum
import numpy as np
import typing

class EnumMeta(enum.EnumMeta):
    def __getitem__(cls, key):
        #if this is an iterable: apply to each value, and construct a tuple
        if isinstance(key, typing.Iterable) and not isinstance(key, str):
            return tuple(map(cls.__getitem__, key))
        #if this is from a flavor scheme: check that it's ours
        if isinstance(key, FlavorScheme):
            if not isinstance(key, cls):
                raise TypeError(f'Value {repr(key)} is not from {cls.__name__} sheme!')
            return key
        #if this is a string find it by name
        if isinstance(key, str):
            try:
                return super().__getitem__(key)
            except KeyError as e:
                    raise KeyError(
                        f'Cannot find key "{key}" in {cls.__name__} sheme! Valid options are {list(cls)}'
                    )
        
        #if this is an int value: find a matching

class FlavorScheme(enum.IntEnum, metaclass=EnumMeta):
    def to_tex(self):
        """LaTeX-compatible string representations of flavor."""
        base = r'\nu'
        if self.is_antineutrino:
            base = r'\overline{\nu}'
        lepton = self.lepton.lower()
        if self.is_muon or self.is_tauon:
            lepton = '\\'+lepton
        return f"${base}_{{{lepton}}}$"
        
    @property
    def is_neutrino(self):
        return not self.is_antineutrino

    @property
    def is_antineutrino(self):
        return '_BAR' in self.name

    @property
    def is_electron(self):
        return self.lepton=='E'

    @property
    def is_muon(self):
        return self.lepton=='MU'

    @property
    def is_tauon(self):
        return self.lepton=='TAU'

    @property
    def is_sterile(self):
        return self.lepton=='S'
        
    @property
    def lepton(self):
        return self.name.split('_')[1]
    

def makeFlavorScheme(name:str, leptons:list):
    enum_class =  FlavorScheme(name, start=0,
                   names = [f'NU_{L}{BAR}' for L in leptons for BAR in ['','_BAR']]
                  )
    return enum_class

TwoFlavor = makeFlavorScheme('TwoFlavor',['E','X'])
ThreeFlavor = makeFlavorScheme('ThreeFlavor',['E','MU','TAU'])
FourFlavor = makeFlavorScheme('ThreeFlavor',['E','MU','TAU','S'])


class FlavorMatrix:
    def __init__(self, 
                 array:np.ndarray,
                 flavors:FlavorScheme,
                 flavors2:FlavorScheme = None
                ):
                    self.array = np.asarray(array)
                    self.flavors1 = flavors
                    self.flavors2 = flavors2 or flavors
                    expected_shape = (len(self.flavors2), len(self.flavors1))
                    if(self.array.shape != expected_shape):
                        raise ValueError(f"FlavorMatrix array shape {self.array.shape} mismatch expected {expected_shape}")

    def _convert_index(self, index):
        if isinstance(index, str) or (not isinstance(index,typing.Iterable)):
            index = [index]
        new_idx = [flavors[idx] for idx,flavors in zip(index, [self.flavors2,self.flavors1])]
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
            data = np.tensordot(self.array, other.array, axes=[1,0])
            return FlavorMatrix(data, other.flavors1, self.flavors2)
        raise TypeError(f"Cannot multiply object of {self.__class__} by {other.__class__}")
    #properties
    @property
    def shape(self):
        return self.array.shape
    @property
    def flavors_right(self):
        return self.flavors1
    @property
    def flavors_left(self):
        return self.flavors2
    @property
    def flavors_from(self):
        return self.flavors_right
    @property
    def flavors_to(self):
        return self.flavors_left
    #methods for creating matrix
    @classmethod
    def zeros(cls, flavors1:FlavorScheme, flavors2:FlavorScheme = None):
        flavors2 = flavors2 or flavors1
        shape = (len(flavors2), len(flavors1))
        data = np.zeros(shape)
        return cls(data, flavors1, flavors2)
        
    @classmethod
    def eye(cls, flavors1:FlavorScheme, flavors2:FlavorScheme = None):
        flavors2 = flavors2 or flavors1
        shape = (len(flavors2), len(flavors1))
        data = np.eye(*shape)
        return cls(data, flavors1, flavors2)

    @classmethod
    def from_function(cls, flavors1:FlavorScheme, flavors2:FlavorScheme = None, *, function):
        flavors2 = flavors2 or flavors1
        data = [[function(f1,f2)
                 for f1 in flavors1]
                for f2 in flavors2]
        return cls(np.array(data,dtype=float), flavors1, flavors2)
    @classmethod
    def identity(cls, flavors1:FlavorScheme, flavors2:FlavorScheme = None):
        return cls.from_function(flavors1, flavors2, lambda f1,f2: 1.*(f1.name==f2.name))


def flavor_matrix(flavor1:FlavorScheme, flavor2:FlavorScheme=None):
    """A decorator for creating the flavor matrix from the given function"""
    flavor2 = flavor2 or flavor1
    def _decorator(func):
        return FlavorMatrix.from_function(flavor1, flavor2, function=func)
    return _decorator

#define the conversion matrices
M_convert = {}

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

#other conversion matrices can be defined