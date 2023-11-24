"""
This module contains utilities to define a new registry supernova model.

Examples usage can be seen in :mod:`snewpy.models.ccsn`.


.. autoclass:: Parameter
    :members:
    
.. autodecorator:: RegistryModel
"""
from functools import wraps
import itertools as it
import inspect
import os
from textwrap import dedent
from warnings import warn
from snewpy._model_downloader import RemoteFileLoader

def deprecated(*names, message='Agrument `{name}` is deprecated'):
    """A function decorator to issue a deprecation warning if a given argument is provided in the wrapped function call.
    
    Parameters
    ----------
    
    names: list of str
        argument names which are deprecated
    message: str
        a template with {name} parameter, to make the message for each argument.

    .. Example::

    @deprecated('foo','bar', message='Argument `{name}` is deprecated and will be removed in SNEWPYv2.0!')
    def test_function(*, foo=1, bar=2, baz=3):
        pass
        
    #calling test_function(foo=1, baz=3) will issue a deprecation warning:
    #    DeprecationWarning: Argument `foo` is deprecated and will be removed in SNEWPYv2.0!
    """
    def _f(func):
        #get function signature
        S = inspect.signature(func)
        @wraps(func)
        def _wrapper(*args, **kwargs):
            #bind signature to find all the parameters
            params = S.bind(*args,**kwargs)
            for name in names:
                if name in params.arguments:
                    warn(message.format(name=name), category=DeprecationWarning, stacklevel=2)
            return func(*args,**kwargs)
        return _wrapper
    return _f


_default_labels={'eos':'EOS', 'magnetic_field_exponent':'B_0 Exponent'}
_default_descriptions={'eos':'Equation of state',
                       'progenitor_mass':'Mass of model progenitor in units Msun',
                       'revival_time':'Time of shock revival in model in units ms',
                       'metallicity':'Progenitor metallicity',
                       'turbmixing_param':'Turbulent mixing parameter alpha_lambda',
                       'rotational_velocity': 'Rotational velocity of progenitor',
                       'magnetic_field_exponent':'Exponent of magnetic field'
                      }
class Parameter(list):
    """A class to describe the model parameter: it's range of allowed values, name, description etc. """
    def __init__(self, values, *,
                 name:str='parameter', 
                 label:str=None, 
                 description:str=None, 
                 desc_values:str=None):
        """
        Parameters
        ----------
        values : iterable
            List of allowed parameter values. If len(values)==1 the parameter is considered fixed (this can be used to add default metadata to the model).
            
        Other Parameters
        ------------------
        name : str
            a variable name for the parameter (used to set the ``label`` and ``description`` if they are `None`)
        label : str
            A key to be used for this parameter in the metadata dictionary (e.g. 'EOS' for name='eos')
            If `None`(default), label will be derived from name
        description : str
            A long description of the parameter, to be used in the docstring generation
            If `None`(default), description will be derived from name
        desc_values:str
            A string representation of the given values range, to be used in the docstring generation
            If `None`(default), desc_values will be just `str(values)`

            
        If label or description is `None` they are derived from the name, 
        by capitalizing and replacing underscores with spaces.
        
        Examples
        --------
        
        >>> # create a parameter object:
        >>> Parameter(range(0,100,1),name='supernova_parameter', desc_values='[0,1,...99]')
        Parameter(name="supernova_parameter", label="Supernova parameter", description="Supernova parameter", values='[0,1,...99]')
        """
        super().__init__(values)
        self.name = name
        self.label = label or _default_labels.get(name, name.replace('_',' ').capitalize())
        self.description = description or _default_descriptions.get(name,self.label)
        self.desc_values = desc_values or str(values)
        
    def __repr__(self):
        return f"{self.__class__.__name__}(name=\"{self.name}\", label=\"{self.label}\", description=\"{self.description}\", values='{self.desc_values}')"
    @property
    def fixed(self):
        """True if this parameter has only one option"""
        return len(self)==1


class ParameterSet:
    """A class to describe all possible combinations of parameter values"""
    def __init__(self, param_validator=None, **params):
        """
        Parameters
        ----------
        params:dict
            Keyword arguments, passing each parameter as iterable or instance of :class:`Parameter`
        param_validator:callable or None
            A function of user parameters (dict), returning true if the passed user parameters are valid.
            If ``None`` (default) - all the combinations are allowed
        """
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

    def fill_default_parameters(self, **user_params)->dict:
        """Set the parameter values if they are missing in 'user_params' and fixed (i.e. have a single option"""
        for name,p in self.params.items():
            if(p.fixed):
                user_params.setdefault(name,p[0])
        return user_params
                
    def validate(self, **user_params):
        """Check that provided user parameters are valid, i.e. are within the given list values"""
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

    def generate_docstring(self, func, **type_annotations)->str:
        #generate docstring
        s = []
        S = inspect.signature(func)
        params = [name for name in self.params.keys() if name in S.parameters]
        for name in params:
            p = self.params[name]
            p_type = type_annotations.get(name,None)
            type_name = p_type.__name__ if p_type else ''
            s+=[f'{name}: {type_name}\n    {p.description}. Valid values are: {p.desc_values}.']
        return '\n'.join(s)
    
def RegistryModel(_init_from_filename=True, _param_validator=None, _config_path='{module_name}.{class_name}',**params):
    """A class decorator for defining the supernova model, that initializes from physics parameters.
    
    Parameters
    -----------
    params:dict
        Keyword arguments, passing each parameter as iterable or instance of :class:`Parameter`
    _param_validator:callable or None
        A function of user parameters (dict), returning true if the passed user parameters are valid.
        If `None` (default) - all the combinations are allowed
    _init_from_filename:bool
        If `True`(default), adds the deprecated initialization from the filename (for the backward compatibility)
    _config_path:str
        config_path : dot-separated path template of the model in the YAML configuration (e.g. "ccsn.Bollig_2016"), or a template.
        Default value is "{module_name}.{class_name}" template. The template parameters will be substituted by:
           "{module_name}" -> the name of the the module where the decorated class is defined
           "{class_name}"  -> the name of the decorated class.
        Ignored, if decorated class has `_config_path` member.
    
    Returns
    -------
        Model class
            A class, inherited from the input class, which has modified `__init__` function (see details below)
    
    Note
    ----
    This decorator helps to construct the daughter of the given class, which automatically implements the constructor with:
    
        * Validation of the input user parameters (see :meth:`RegistryModel.validate`)
        * Generated constructor docstring based on allowed parameters
        * Populates the `self.metadata` from the given user parameters
        * Optional (deprecated) "filename" argument, and calls initialization from the filename

        The decorated class:
            * *must* inherit from the loader class
            * *must* implement the __init__ method, which
            
                1. defines the filename from the given user parameters
                2. optionally modify the `self.metadata` dictionary, and 
                3. calls the corresponding loader class constructor
            * *can* define the ``_metadata_from_filename (self, filename:str)->dict``, to help populate the metadata when a `filename` argument is passed
    
    """
    pset:ParameterSet = ParameterSet(param_validator=_param_validator, **params)
    def _wrap(base_class):
        class c(base_class, RemoteFileLoader):
            _config_path="{module_path}.{class_name}"
            _doc_params_ = {
                'Other parameters': pset.generate_docstring(base_class.__init__, **base_class.__init__.__annotations__),
                'Raises':"""
                FileNotFoundError
                    If a file for the chosen model parameters cannot be found
                ValueError:
                    If a combination of parameters is invalid when loading from parameters"""
            }
            @classmethod
            def _generate_docstring(cls)->str:
                docstring = dedent(base_class.__init__.__doc__ or '')
                for section, desc in cls._doc_params_.items():
                    s = f'{section}\n'+'-'*len(section)+'\n'+dedent(desc)
                    docstring+='\n'+s+'\n'
                return docstring
                
            def __init__(self, **kwargs):
                # enforce the default parameters
                S = inspect.signature(self.__init__)
                params = S.bind(**kwargs)
                params.apply_defaults()
                #print(f"Calling init with params={params}")
                kwargs = params.kwargs
                #select the parameters which correspond to metadata
                param_kwargs = {name:val for name,val in kwargs.items() if name in pset.params}
                param_kwargs = pset.fill_default_parameters(**param_kwargs)
                # validate the input parameters
                pset.validate(**param_kwargs)
                #Store model metadata
                self.metadata = {pset[name].label: value for name,value in param_kwargs.items()}
                #call the constructor with only the needed arguments (the rest are only for metadata)
                S = inspect.signature(super().__init__)
                init_params = {name:val for name,val in kwargs.items() if name in S.parameters}
                print(init_params)
                return super().__init__(**init_params)
                
            @classmethod
            def get_param_combinations(cls)->tuple:
                """Get all valid combinations of parameters for a this registry model.
                Returns
                -------
                valid_combinations: tuple[dict]
                    A tuple of all valid parameter combinations stored as Dictionaries"""
                return pset.valid_combinations_dict
        #generate the configuration path
        if not hasattr(c,"_config_path"):
            c._config_path="{module_name}.{class_name}"
        #generate the configuration path
        c._config_path=_config_path.format(module_name=base_class.__module__.split('.')[-1], class_name=base_class.__name__)
        print(f'config_path = {c._config_path}')
        #generate the docstring
        c.__doc__ = base_class.__doc__
        c.__init__.__doc__ = c._generate_docstring()
        #copy the signature
        c.__init__.__signature__ = inspect.signature(base_class.__init__)
        if not _init_from_filename:
            c.__qualname__ = base_class.__qualname__
            c.__name__ = base_class.__name__
            c.__module__ = base_class.__module__
            return c

        class c1(c):
            @deprecated('filename')
            def __init__(self, filename:str=None, **kwargs):
                if filename is not None:
                    self.metadata = {}
                    if hasattr(self,'_metadata_from_filename'):
                        self.metadata = self._metadata_from_filename(filename)
                    super(base_class, self).__init__(filename=os.path.abspath(filename), metadata=self.metadata)
                else:
                    super().__init__(**kwargs)
        #generate the docstring
        c1.__doc__ = c.__doc__
        c1._doc_params_ = {'Parameters':"""filename: str\n    Absolute or relative path to the file with model data. This argument is deprecated.""", **c._doc_params_}
        c1.__init__.__doc__ = c1._generate_docstring()
        
        #update the call signature
        S = inspect.signature(c)
        S1 = inspect.signature(c1.__init__)
        #set default values to None if they are not set
        kw_params = [p.replace(kind=inspect.Parameter.KEYWORD_ONLY) for name,p in S.parameters.items()]
        params = [S1.parameters['self'],S1.parameters['filename'],*kw_params]
        c1.__init__.__signature__ = S.replace(parameters=params)
        c1.__qualname__ = base_class.__qualname__
        c1.__name__ = base_class.__name__
        c1.__module__ = base_class.__module__
        return c1
        
    return _wrap    