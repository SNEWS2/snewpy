from functools import wraps
import itertools as it
import inspect
import os

_default_labels={'eos':'EOS'}
_default_descriptions={'eos':'Equation of state',
                       'progenitor_mass':'Mass of model progenitor in units Msun',
                       'revival_time':'Time of shock revival in model in units ms',
                       'metallicity':'Progenitor metallicity',
                       'turbmixing_param':'Turbulent mixing parameter alpha_lambda',
                       'rotational_velocity': 'Rotational velocity of progenitor',
                       'magnetic_field_exponent':'Exponent of magnetic field'
                      }
class Parameter(list):
    def __init__(self, values, *, name:str='parameter', label=None, description=None, desc_values=None):
        super().__init__(values)
        self.name = name
        self.label = label or _default_labels.get(name, name.replace('_',' ').capitalize())
        self.description = description or _default_descriptions.get(name,self.label)
        self.desc_values = desc_values or str(values)
        
    def __repr__(self):
        return f"{self.__class__.__name__}(name=\"{self.name}\", label=\"{self.label}\", description=\"{self.description}\", values='{self.desc_values}')"


class ParameterSet:
    def __init__(self, param_validator=None, **params):
        #construct Parameter classes from given arrays
        self.params = {}
        for name,val in params.items():
            if not isinstance(val,Parameter):
                val = Parameter(values=val,name=name)
            self.params[name]=val
                
        self.combinations = list(it.product(*[list(v) for v in self.params.values()]))
        self.valid_combinations = []
        if param_validator:
            for p in self.combinations:
                pdict = dict(zip(params.keys(),p))
                is_valid = param_validator(pdict)
                if is_valid:
                    self.valid_combinations+=[p]
        else:
            self.valid_combinations = self.combinations
        self.valid_combinations_dict = [dict(zip(self.params,p)) for p in self.valid_combinations]
        
    def validate(self, **user_params):
        #check that we have all correct parameters
        for name in user_params:
            if name not in self.params:
                raise ValueError(f"Unexpected parameter '{name}', allowed parameters are {list(self.params.keys())}")
            if user_params[name] not in self.params[name]:
                #check if all the parameters are in allowed ranges
                raise ValueError(f"Invalid parameter value for '{name}'={user_params[name]}. Allowed values are: {self.params[name].desc_values}")
        for name in self.params:
            if name not in user_params:
                raise ValueError(f"Missing parameter '{name}'")
        ptuple = tuple([user_params[name] for name in self.params])
        if not ptuple in self.valid_combinations:
            valid_dicts = [dict(zip(self.params.keys(),pars)) for pars in self.valid_combinations]
            raise ValueError(f"Invalid parameters combination: {user_params}")

    def __getitem__(self, name:str):
        return self.params.__getitem__(name)

    def __repr__(self):
        s = f"{self.__class__.__name__}:\n"
        s+='\n'.join([f' * \t{name}={ps}' for name,ps in self.params.items()])
        return s

    def generate_docstring(self, **type_annotations)->str:
        #generate docstring
        s = []
        for name,p in self.params.items():
            p_type = type_annotations.get(name,None)
            type_name = p_type.__name__ if p_type else ''
            s+=[f'{name}: {type_name}\n    {p.description}. Valid values are: {p.desc_values}.']
        return '\n'.join(s)


def RegistryModel(_init_from_filename=True, _param_validator=None, **params):
    pset:ParameterSet = ParameterSet(param_validator=_param_validator, **params)
    def _wrap(base_class):
        class c(base_class):
            def __init__(self, **kwargs):
                """
Keyword parameters
------------------
{PARAMETERS} 

Raises
------
ValueError
    If a combination of parameters is invalid when loading from parameters"""
                # validate the input parameters                  
                pset.validate(**kwargs)
                # Store model metadata.
                self.metadata = {pset[name].label: value for name,value in kwargs.items()}
                return super().__init__(**kwargs)
                
            @classmethod
            def get_param_combinations(cls)->tuple:
                """Get all valid combinations of parameters for a this registry model.
                Returns
                -------
                valid_combinations: tuple[dict]
                    A tuple of all valid parameter combinations stored as Dictionaries"""
                return pset.valid_combinations_dict
        c.__init__.__doc__ = c.__init__.__doc__.format(
            PARAMETERS=pset.generate_docstring(**base_class.__init__.__annotations__)
        )
        c.__init__.__signature__ = inspect.signature(base_class.__init__)
        if not _init_from_filename:
            c.__name__ = base_class.__name__
            return c

        class c1(c):
            def __init__(self, filename:str=None, **kwargs):
                """
Parameters
----------
filename: str
    Absolute or relative path to the file with model data. This argument is deprecated.
    """
                #deprecated initialization from the 
                print(base_class)
                if filename is not None:
                    self.metadata = {}
                    if hasattr(self,'_metadata_from_filename'):
                        self.metadata = self._metadata_from_filename(filename)
                    super(base_class, self).__init__(filename=os.path.abspath(filename), metadata=self.metadata)
                else:
                    super().__init__(**kwargs)
                    
        c1.__init__.__doc__+=c.__init__.__doc__   
        #update the call signature
        S = inspect.signature(c)
        S1 = inspect.signature(c1.__init__)
        params = [S1.parameters['self'],S1.parameters['filename'],
                  *(p.replace(kind=inspect.Parameter.KEYWORD_ONLY, default=None) for name,p in S.parameters.items())
                 ]
        c1.__init__.__signature__ = S.replace(parameters=params)
        c1.__name__ = base_class.__name__
        return c1
        
    return _wrap    