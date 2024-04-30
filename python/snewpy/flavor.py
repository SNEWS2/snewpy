def MakeFlavorScheme(name:str, leptons:list[str]):
    return enum.Enum(name, start=0, type=FlavorScheme,
                   names = [f'NU_{L}{BAR}' for L in leptons for BAR in ['','_BAR']]
                  )

TwoFlavor = MakeFlavorScheme('TwoFlavor',['E','X'])
ThreeFlavor = MakeFlavorScheme('ThreeFlavor',['E','MU','TAU'])
FourFlavor = MakeFlavorScheme('ThreeFlavor',['E','MU','TAU','S'])

class FlavorScheme:
    def to_tex(self):
        """LaTeX-compatible string representations of flavor."""
        if self.is_antineutrino:
            return r'$\overline{{\nu}}_{0}$'.format(self.name[3].lower())
        return r'$\{0}$'.format(self.name.lower())

    @classmethod
    def to_value(cls,key)->int:
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
    return enum.Enum(name, start=0, type=FlavorScheme,
                   names = [f'NU_{L}{BAR}' for L in leptons for BAR in ['','_BAR']]
                  )

TwoFlavor = MakeFlavorScheme('TwoFlavor',['E','X'])
ThreeFlavor = MakeFlavorScheme('ThreeFlavor',['E','MU','TAU'])
FourFlavor = MakeFlavorScheme('ThreeFlavor',['E','MU','TAU','S'])


class FlavorMatrix:
    def __init__(self, 
                 array:np.ndarray,
                 flavors:FlavorScheme,
                 flavors2:FlavorScheme|None = None
                ):
                    self.array = np.asarray(array)
                    self.flavors1 = flavors
                    self.flavors2 = flavors2 or flavors
                    assert self.array.shape == (len(self.flavors2), len(self.flavors1))

    def _convert_index(self, index):
        if isinstance(index, str) or (not isinstance(index,typing.Iterable)):
            index = [index]
        new_idx = [flavors.to_value(idx) for idx,flavors in zip(index, [self.flavors2,self.flavors1])]
        print(new_idx)
        return tuple(new_idx)
        
    def __getitem__(self, index):
        return self.array[self._convert_index(index)]
        
    def __setitem__(self, index, value):
        self.array[self._convert_index(index)]=value
        
    def __repr__(self):
        s = f'{self.__class__.__name__}:<{self.flavors1.__name__}->{self.flavors2.__name__}>:'
        s+='\n'+repr(self.array)
        return s
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
    def identity(cls, flavors1:FlavorScheme, flavors2:FlavorScheme|None = None):
        flavors2 = flavors2 or flavors1
        data = [[f1.name==f2.name 
                 for f1 in flavors1]
                for f2 in flavors2]
        return cls(np.array(data,dtype=float), flavors1, flavors2)