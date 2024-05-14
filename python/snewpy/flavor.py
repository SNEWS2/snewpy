import enum
import numpy as np
import typing

class EnumMeta(enum.EnumMeta):
    def __getitem__(cls, key):
        #if this is an iterable: apply to each value, and construct a tuple
        if isinstance(key, slice):
            return slice(cls[key.start],cls[key.stop],key.step)
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
        
        if key is None:
            return None
        #if this is anything else - treat it as a slice
        return np.array(list(cls.__members__.values()),dtype=object)[key]

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

    @classmethod
    def from_lepton_names(cls, name:str, leptons:list):
        enum_class =  cls(name, start=0, names = [f'NU_{L}{BAR}' for L in leptons for BAR in ['','_BAR']])
        return enum_class
        
    @classmethod
    def take(cls, index):
        return cls[index]

TwoFlavor = FlavorScheme.from_lepton_names('TwoFlavor',['E','X'])
ThreeFlavor = FlavorScheme.from_lepton_names('ThreeFlavor',['E','MU','TAU'])
FourFlavor = FlavorScheme.from_lepton_names('FourFlavor',['E','MU','TAU','S'])

class FlavorMatrix:
    def __init__(self, 
                 array:np.ndarray,
                 flavor:FlavorScheme,
                 from_flavor:FlavorScheme = None
                ):
                    self.array = np.asarray(array)
                    self.flavor_out = flavor
                    self.flavor_in = from_flavor or flavor
                    expected_shape = (len(self.flavor_out), len(self.flavor_in))
                    if(self.array.shape != expected_shape):
                        raise ValueError(f"FlavorMatrix array shape {self.array.shape} mismatch expected {expected_shape}")
       
    def _convert_index(self, index):
        if isinstance(index, str) or (not isinstance(index,typing.Iterable)):
            index = [index]
        new_idx = [flavors[idx] for idx,flavors in zip(index, [self.flavor_out, self.flavor_in])]
        return tuple(new_idx)
        
    def __getitem__(self, index):
        return self.array[self._convert_index(index)]
        
    def __setitem__(self, index, value):
        self.array[self._convert_index(index)] = value

    def _repr_short(self):
        return f'{self.__class__.__name__}:<{self.flavor_in.__name__}->{self.flavor_out.__name__}> shape={self.shape}'
        
    def __repr__(self):
        s=self._repr_short()+'\n'+repr(self.array)
        return s
    def __eq__(self,other):
        return self.flavor_in==other.flavor_in and self.flavor_out==other.flavor_out and np.allclose(self.array,other.array)
                
    def __matmul__(self, other):
        if isinstance(other, FlavorMatrix):
            try:
                data = np.tensordot(self.array, other.array, axes=[1,0])
                return FlavorMatrix(data, self.flavor_out, from_flavor = other.flavor_in)
            except Exception as e:
                raise ValueError(f"Cannot multiply {self._repr_short()} by {other._repr_short()}") from e
        elif hasattr(other, '__rmatmul__'):
            return other.__rmatmul__(self)        
        raise TypeError(f"Cannot multiply object of {self.__class__} by {other.__class__}")
    #properties
    @property
    def shape(self):
        return self.array.shape
    @property
    def flavor(self):
        return self.flavor_out
        
    @classmethod
    def eye(cls, flavor:FlavorScheme, from_flavor:FlavorScheme = None):
        from_flavor = from_flavor or flavor
        shape = (len(from_flavor), len(flavor))
        data = np.eye(*shape)
        return cls(data, flavor, from_flavor)

    @classmethod
    def from_function(cls, flavor:FlavorScheme, from_flavor:FlavorScheme = None):
        """A decorator for creating the flavor matrix from the given function"""
        from_flavor = from_flavor or flavor
        def _decorator(function):
            data = [[function(f1,f2)
                     for f2 in from_flavor]
                    for f1 in flavor]
                    
            return cls(np.array(data,dtype=float),  flavor, from_flavor)
        return _decorator
    #flavor conversion utils
    
def conversion_matrix(from_flavor:FlavorScheme, to_flavor:FlavorScheme):
    @FlavorMatrix.from_function(to_flavor, from_flavor)
    def convert(f1, f2):
        if f1.name == f2.name:
            return 1.
        if (f1.is_neutrino == f2.is_neutrino) and (f2.lepton == 'X' and f1.lepton in ['MU', 'TAU']):
            # convert from TwoFlavor to more flavors
            return 1.
        if (f1.is_neutrino == f2.is_neutrino) and (f1.lepton == 'X' and f2.lepton in ['MU', 'TAU']):
            # convert from more flavors to TwoFlavor
            return 0.5
        return 0.
    return convert

FlavorScheme.conversion_matrix = classmethod(conversion_matrix)
EnumMeta.__rshift__ = conversion_matrix
EnumMeta.__lshift__ = lambda f1,f2:conversion_matrix(f2,f1)
