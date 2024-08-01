import enum
import numpy as np
import typing
import snewpy.utils

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
            #prepare the string
            key=key.upper()
            #check cases for neutrinos and antineutrinos
            if key=='NU':
                return tuple([f for f in cls.__members__.values() if f.is_neutrino])
            elif key=='NU_BAR':
                return tuple([f for f in cls.__members__.values() if not f.is_neutrino])
            #add the prefix if needed
            if not key.startswith('NU_'):
                key = 'NU_'+key
            #try to get the proper values
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
                    if(self.array.shape[:2] != expected_shape):
                        raise ValueError(f"FlavorMatrix array shape {self.array.shape} mismatch expected {expected_shape}")
       
    def _convert_index(self, index):
        if isinstance(index, str) or (not isinstance(index,typing.Iterable)):
            index = [index]
        #convert flavor dimensions
        new_idx = [flavors[idx] for idx,flavors in zip(index, [self.flavor_out, self.flavor_in])]
        #try to convert first index to list of bracketed index
        try:
            new_idx[0] = [[f] for f in new_idx[0]]
        except:
            #it was not iterable
            pass 
        #add remaining dimensions
        new_idx+=list(index[2:])
        return tuple(new_idx)
        
    def __getitem__(self, index):
        index = self._convert_index(index)
        array = self.array[index]
        #if requested all the flavors, then return a FlavorMatrix using the remaining index
        if index[:2]==(slice(None),slice(None)):
            return FlavorMatrix(array, self.flavor_out, self.flavor_in)
        else:
            return array
        
    def __setitem__(self, index, value):
        self.array[self._convert_index(index)] = value

    def _repr_short(self):
        return f'{self.__class__.__name__}:<{self.flavor_in.__name__}->{self.flavor_out.__name__}> shape={self.shape}'
        
    def __repr__(self):
        s=self._repr_short()
        if(len(self.shape)==2):
            s+='\n'+repr(self.array)
        return s
        
    def _repr_markdown_(self):
        if(len(self.shape)>2):
            return self.__repr__()
        #make a markdown table
        res = [f'**{self.__class__.__name__}**:<`{self.flavor_in.__name__}`->`{self.flavor_out.__name__}`> shape={self.shape}','']
        flavors0 = [f.to_tex() for f in self.flavor_in]
        flavors1 = [f.to_tex() for f in self.flavor_out]
        res+=['|'.join(['*']+flavors1)]
        res+=['|'.join(['-:',]*(len(flavors1)+1))]
        for f0 in self.flavor_in:
            line = [f0.to_tex()]+[f'{v:.3g}' for v in self[:,f0]]
            res+=['|'.join(line)]
        return '\n'.join(res)
    
    def __eq__(self,other):
        return self.flavor_in==other.flavor_in and self.flavor_out==other.flavor_out and np.allclose(self.array,other.array)
    @property
    def real(self):
        return FlavorMatrix(self.array.real, self.flavor_out, from_flavor = self.flavor_in)
    @property
    def imag(self):
        return FlavorMatrix(self.array.imag, self.flavor_out, from_flavor = self.flavor_in)
    
    def abs(self):
        return FlavorMatrix(np.abs(self.array), self.flavor_out, from_flavor = self.flavor_in)
    def abs2(self):
        return FlavorMatrix(np.abs(self.array**2), self.flavor_out, from_flavor = self.flavor_in)
        
    def __mul__(self, other):
        if isinstance(other, FlavorMatrix):
            if not ((other.flavor_in==self.flavor_in)and(other.flavor_out==self.flavor_out)):
                raise TypeError(f"Cannot multiply matrices with different flavor schemes: {self._repr_short()} and {other._repr_short()}")
            other = other.array
        return FlavorMatrix(self.array*other, self.flavor_out, from_flavor = self.flavor_in)
    def __matmul__(self, other):
        if isinstance(other, FlavorMatrix):
            assert self.flavor_in==other.flavor_out, f"Incompatible spaces {self.flavor_in}!={other.flavor_out}"
            try:
                m0, m1 = self.array, other.array
                ndims = max(m0.ndim, m1.ndim)
                m0,m1 = [snewpy.utils.expand_dimensions_to(m, ndim=ndims) for m in [m0,m1]]
                array = np.einsum('ij...,jk...->ik...',m0,m1)
                np.tensordot(self.array, other.array, axes=[1,0])
                return FlavorMatrix(array, self.flavor_out, from_flavor = other.flavor_in)
            except Exception as e:
                raise ValueError(f"Cannot multiply {self._repr_short()} by {other._repr_short()}") from e
        elif hasattr(other, '__rmatmul__'):
            return other.__rmatmul__(self)
        elif isinstance(other, dict):
            #try to multiply to the dict[Flavor:flux]
            #we assume that it has the same Flavor scheme!
            array =  np.stack([other[flv] for flv in self.flavor_in])
            result = np.tensordot(self.array, array, axes=[1,0])
            return {flv:result[n] for n,flv in enumerate(self.flavor_out)}
        raise TypeError(f"Cannot multiply object of {self.__class__} by {other.__class__}")
    #properties
    @property
    def shape(self):
        return self.array.shape
    @property
    def flavor(self):
        return self.flavor_out
    @property
    def T(self):
        return self.transpose()
    
    def transpose(self):
        "transposed version of the matrix: reverse the flavor dimensions"
        return FlavorMatrix(array = self.array.swapaxes(0,1), 
                            flavor = self.flavor_in, from_flavor=self.flavor_out)
    def conjugate(self):
        "apply complex conjugate"
        return FlavorMatrix(array = self.array.conjugate(), 
                            flavor = self.flavor_out, from_flavor=self.flavor_in)

    @classmethod
    def eye(cls, flavor:FlavorScheme, from_flavor:FlavorScheme = None, extra_dims=[]):
        return cls.from_function(flavor,from_flavor)(lambda f1,f2: f1==f2*np.ones(shape=extra_dims))
    @classmethod
    def zeros(cls, flavor:FlavorScheme, from_flavor:FlavorScheme = None, extra_dims=[]):
        from_flavor = from_flavor or flavor
        shape = (len(from_flavor), len(flavor), *extra_dims)
        data = np.zeros(shape)
        return cls(data, flavor, from_flavor)
    @classmethod
    def from_function(cls, flavor:FlavorScheme, from_flavor:FlavorScheme = None):
        """A decorator for creating the flavor matrix from the given function"""
        if from_flavor is None: 
            from_flavor = flavor
        def _decorator(function):
            data = [[function(f1,f2)
                     for f2 in from_flavor]
                    for f1 in flavor]
                    
            return cls(np.array(data,dtype=float),  flavor, from_flavor)
        return _decorator
    #flavor conversion utils
    def convert_to_flavor(self, flavor_out:FlavorScheme|None=None, flavor_in:FlavorScheme|None=None):
        if flavor_out is None and flavor_in is None:
            raise ArgumentError('Provide flavor_in and/or flavor_out argument')
        M = self
        if flavor_in is not None:
            M = M@(M.flavor_in<<flavor_in)
        if flavor_out is not None:
            M = (flavor_out<<M.flavor_out)@M
        return M
    
    def __rshift__(self, flavor:FlavorScheme):
        return self.convert_to_flavor(flavor_out=flavor)
    def __lshift__(self, flavor:FlavorScheme):
        return self.convert_to_flavor(flavor_in=flavor)
        
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
    #check if the conversion matrix is not zero
    if np.allclose(convert.array, 0):
        raise RuntimeError(f"Conversion matrix {convert._repr_short()} is zero!")
    return convert

def rshift(flv:FlavorScheme, obj:FlavorScheme|FlavorMatrix)->FlavorMatrix:
    if isinstance(obj, EnumMeta):
        return conversion_matrix(from_flavor=flv,to_flavor=obj)
    elif hasattr(obj, '__lshift__'):
        return obj<<flv
    else:
        raise TypeError(f'Cannot apply flavor conversion to object of type {type(obj)}')

def lshift(flv:FlavorScheme, obj:FlavorScheme|FlavorMatrix)->FlavorMatrix:
    if isinstance(obj, EnumMeta):
        return conversion_matrix(from_flavor=obj,to_flavor=flv)
    elif hasattr(obj, '__rshift__'):
        return obj>>flv
    else:
        raise TypeError(f'Cannot apply flavor conversion to object of type {type(obj)}')
        
        
FlavorScheme.conversion_matrix = classmethod(conversion_matrix)
EnumMeta.__rshift__ = rshift
EnumMeta.__lshift__ = lshift
