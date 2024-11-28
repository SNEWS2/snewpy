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
import numpy as np
import re
from snewpy._model_downloader import RegistryFileLoader

from astropy import table
from tqdm.auto import tqdm

all_models = set()
    
def get_models_table(init:bool=False)->table.QTable:
    """Get the astropy.Table with all the possible Metadata parameters.

    Parameters
    ----------
    init:bool
        If true - initialize each model to obtain the runtime metadata parameters (takes longer)
        Default is false, i.e. use only input parameters
    """
    #prepare the table
    tables = {}
    for model in tqdm(all_models):
        param_combinations = model.get_param_combinations()
        if(init):
            metadata = [model(**params).metadata for params in param_combinations]
        else:
            metadata = [{model.parameters[key].label:val for key,val in params.items()} for params in param_combinations]
            
        md_table = table.QTable(metadata)
        md_table['model']='.'.join([model.__module__,model.__name__])
        md_table['init_params']=param_combinations
        #set the fixed metadata
        if hasattr(model,"__metadata__"):
            for key,value in model.__metadata__.items():
                md_table[key]=value
        tables[model]=md_table
    
    df = table.vstack(list(tables.values()))
    return df
    
def _can_decorate_class_or_func(func_decorator):
    """make the function decorator act as class decorator:
    if decorated object is a class, wrap its "__init__" function.
    """
    def _wrapper(obj):
        if inspect.isclass(obj):
            obj.__init__=func_decorator(obj.__init__)
            return obj
        else:
            return func_decorator(obj)
    return _wrapper
            
def deprecated(*names, message='Argument `{name}` is deprecated.'):
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
    @_can_decorate_class_or_func
    def _wrapper(func):
        #get function signature
        S = inspect.signature(func)
        @wraps(func)
        def _f(*args, **kwargs):
            #bind signature to find all the parameters
            params = S.bind(*args,**kwargs)
            for name in names:
                if name in params.arguments:
                    warn(message.format(name=name), category=UserWarning, stacklevel=2)
            return func(*args,**kwargs)
        return _f
    return _wrapper

def map_arguments(**names_dict):
    """map function arguments to create a function with new signature,
    and with updated docstring
    
    Parameters
    ----------
    names_dict: dict
        Argument names in format {old_name:new_name}
    
    Examples
    --------
    >>> @map_arguments(bar='BAR')
    >>> def foo(bar, baz):
    >>>     print(f'bar={bar} baz={baz}')
    >>>
    >>> foo(BAR=123, baz=10) 
    >>> #bar=123 baz=10
    """
    def _wrapper(func):
        S = inspect.signature(func) #old signature
        old2new = names_dict
        new2old = {v:k for k,v in old2new.items()}
        
        @wraps(func)
        def func_new(*args, **kwargs):
            params = S1.bind(*args,**kwargs)
            kwargs_old = params.kwargs
            kwargs_new = {new2old.get(name,name):val for name,val in kwargs_old.items()}
            return func(*params.args,**kwargs_new)
            
        #update signature
        S1 = S.replace(parameters=[
            p.replace(name=old2new.get(p.name,p.name)) for p in S.parameters.values()
            ])
        func_new.__signature__=S1
        #update the docstring
        def _update_docs(docs):
            for name_old,name_new in old2new.items():
                docs = re.sub(rf'\n{name_old}([\s:].*\n)',rf'\n{name_new}\g<1>',docs)
            return docs
        func_new.__doc__ = _update_docs(func.__doc__)
        if(hasattr(func_new,'__init__')):
            func_new.__init__.__doc__ = _update_docs(func.__init__.__doc__)
        return func_new
    return _wrapper
    
def _expand_defaults(func, **defaults):
    """Update the signature of the given function with the parameters with default values.
    Examples
    --------
    >>> def foo(x, y=1):
    >>>     pass
    >>> foo?
    Signature: foo(x, y=1)
    >>> _expand_defaults(foo, x=123, z=456)
    >>> foo?
    Signature: foo(x=123, y=500, *, z=456)
    """
    S = inspect.signature(func)
    params = list(S.parameters.values())
    #defaults.update(**{p.name:p.default for p in params if p.default!=p.empty})
    params_updated = [p.replace(default=defaults.pop(p.name, p.default)) for p in params]
    #add the remaining default parameters
    params_updated+=[inspect.Parameter(name,kind=inspect.Parameter.KEYWORD_ONLY, default=value) for name,value in defaults.items()]
    func.__signature__ = S.replace(parameters=params_updated)

_default_labels={'eos':'EOS', 'magnetic_field_exponent':'B_0 Exponent'}
_default_descriptions={'eos':'Equation of state',
                       'progenitor_mass':'Mass of model progenitor in units Msun',
                       'revival_time':'Time of shock revival in model in units ms',
                       'metallicity':'Progenitor metallicity',
                       'turbmixing_param':'Turbulent mixing parameter alpha_lambda',
                       'rotational_velocity': 'Rotational velocity of progenitor',
                       'magnetic_field_exponent':'Exponent of magnetic field'
                      }
class Parameter:
    """A class to describe the model parameter: it's range of allowed values, name, description etc. """
    def __init__(self, values, *,
                 name:str=None, 
                 label:str=None, 
                 description:str=None, 
                 desc_values:str=None,
                 precision:int=None):
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
        precision:int or None
            Number of decimals, used to round the input value if when testing if given element is in the list of allowed values. This should be used only if the `values` are an array of floats or `Quantity`.
            If `None` (default) no rounding is applied.
            
        If label or description is `None` they are derived from the name, 
        by capitalizing and replacing underscores with spaces.
        
        Examples
        --------
        
        >>> # create a parameter object:
        >>> Parameter(range(0,100,1),name='supernova_parameter', desc_values='[0,1,...99]')
        Parameter(name="supernova_parameter", label="Supernova parameter", description="Supernova parameter", values='[0,1,...99]')
        """
        self.values = values
        self.label = label 
        self.description = description
        self.desc_values = desc_values or str(values)
        self.precision = precision
        if(name):
            self.set_description_and_label(name)
        
    def set_description_and_label(self, name:str):
        """Generate the description and label from the parameter given name"""
        self.label = self.label or _default_labels.get(name, name.replace('_',' ').capitalize())
        self.description = self.description or _default_descriptions.get(name,self.label)
        
    def apply_precision(self, value):
        if self.precision:
            if hasattr(self.values, 'unit'):
                #convert unit before 
                value = value<<self.values.unit 
            value = np.around(value,decimals=self.precision)
        return value
        
    def __contains__(self, value):
        value = self.apply_precision(value)
        return value in self.values
    def __iter__(self):
        return self.values.__iter__()
    def __len__(self):
        return self.values.__len__()
    def __getitem__(self, index):
        return self.values.__getitem__(index)
    def __repr__(self):
        return f"{self.__class__.__name__}(label=\"{self.label}\", description=\"{self.description}\", values={self.desc_values})"
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
            else:
                val.set_description_and_label(name)
            self.params[name]=val
        #fill the valid parameter combinations
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
        #fill the default values for fixed parameters
        self.defaults = {name:p[0] for name,p in self.params.items() if p.fixed}
            
    def _fill_default_parameters(self, **user_params)->dict:
        """Set the parameter values if they are missing in 'user_params' and fixed (i.e. have a single option"""
        return dict(self.defaults, **user_params)

    def _apply_precision(self, **user_params)->dict:
        """Check that we have all required parameters in `user_params`.
        Apply parameter precision and unit conversion for every given parameter"""
        for name,val in user_params.items():
            if name in self.params:
                user_params[name]=self.params[name].apply_precision(val)
        return user_params
        
    def _check_all_parameters_present(self, **user_params):
        """Check that `user_params` contains all needed parameters, and nothing extra"""
        #check that we have all correct parameters
        for name in user_params:
            if name not in self.params:
                raise ValueError(f"Unexpected parameter '{name}', allowed parameters are {list(self.params.keys())}")
        for name in self.params:
            if name not in user_params:
                raise ValueError(f"Missing parameter '{name}'")
                
    def _check_parameter_values(self, **user_params):
        """Check that provided user parameters are valid, i.e. are within the given list values"""
        for name in user_params:
            if user_params[name] not in self.params[name]:
                raise ValueError(f"Invalid parameter value for '{name}'={user_params[name]}. Allowed values are: {self.params[name].desc_values}")
        #check that parameter combination is valid
        ptuple = tuple([user_params[name] for name in self.params])
        if not ptuple in self.valid_combinations:
            valid_dicts = [dict(zip(self.params.keys(),pars)) for pars in self.valid_combinations]
            raise ValueError(f"Invalid parameters combination: {user_params}")
            
    def validate(self, **user_params):
        """Check that provided user parameters are valid, i.e. are within the given list values"""
        #check that we have all correct parameters
        self._check_all_parameters_present(**user_params)
        user_params = self._apply_precision(**user_params)
        self._check_parameter_values(**user_params)

    def __getitem__(self, name:str):
        return self.params.__getitem__(name)
    def items(self):
        return self.params.items()
    def values(self):
        return self.params.values()
    def __iter__(self):
        return self.params.__iter__()
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
    
def RegistryModel(_param_validator=None, **params):
    """A class decorator for defining the supernova model, that initializes from physics parameters.
    
    Parameters
    -----------
    params:dict
        Keyword arguments, passing each parameter as iterable or instance of :class:`Parameter`
    _param_validator:callable or None
        A function of user parameters (dict), returning true if the passed user parameters are valid.
        If `None` (default) - all the combinations are allowed

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

    If a parameter has a single allowed value (i.e. fixed), this value will be default for this parameter in constructor.
    
    """
    pset:ParameterSet = ParameterSet(param_validator=_param_validator, **params)
    def _wrap(base_class):
        class c(base_class, RegistryFileLoader):
            parameters:ParameterSet = pset
            @classmethod
            def _generate_docstring(cls)->str:
                docstring = dedent(base_class.__init__.__doc__ or '')
                for section, desc in cls._doc_params_.items():
                    s = f'{section}\n'+'-'*len(section)+'\n'+dedent(desc)
                    docstring+='\n'+s+'\n'
                return docstring
                
            def __init__(self, *args, **kwargs):
                # enforce the default parameters
                params = inspect.signature(self.__init__).bind(*args, **kwargs)
                params.apply_defaults()
                #select the parameters which correspond to metadata
                arguments = {name:val for name,val in params.arguments.items() if name in self.parameters}
                arguments = self.parameters._fill_default_parameters(**arguments)
                arguments = self.parameters._apply_precision(**arguments)
                # validate the input parameters
                self.parameters.validate(**arguments)
                #Store model metadata
                self.metadata = {self.parameters[name].label: value for name,value in arguments.items()}
                #call the constructor with only the needed arguments (the rest are only for metadata)
                S = inspect.signature(super().__init__)
                init_params = {name:val for name,val in arguments.items() if name in S.parameters}
                return super().__init__(**init_params)
                
            @classmethod
            def get_param_combinations(cls)->tuple:
                """Get all valid combinations of parameters for a this registry model.
                Returns
                -------
                valid_combinations: tuple[dict]
                    A tuple of all valid parameter combinations stored as Dictionaries"""
                return cls.parameters.valid_combinations_dict
            
            param = {name: p.values for name, p in parameters.items()}
        #fill the constructor signature
        c.__init__.__signature__ = inspect.signature(base_class.__init__)
        #If we have "fixed" parameters in ParameterSet, add them to the signature as keyword arguments
        defaults = c.parameters.defaults
        _expand_defaults(c.__init__, **defaults)
        #generate the docstring
        c._doc_params_ = {
                'Other parameters': pset.generate_docstring(c.__init__, **c.__init__.__annotations__),
                'Raises':"""
                FileNotFoundError
                    If a file for the chosen model parameters cannot be found
                ValueError
                    If a combination of parameters is invalid when loading from parameters"""
            }
        c.__doc__ = base_class.__doc__
        c.__init__.__doc__ = c._generate_docstring()
        #fill the class and module name to be the same as in class
        c.__qualname__ = base_class.__qualname__
        c.__name__ = base_class.__name__
        c.__module__ = base_class.__module__
        #set the configuration path
        if not hasattr(c, '_config_path'):
            module_name = c.__module__.split(".")[-1]
            c._config_path = f'{module_name}.{c.__name__}'
        #register the model in the list
        global all_models
        all_models.add(c)
        return c
    return _wrap

def legacy_filename_initialization(c):
    """Wrap the model class, adding a filename argument in the init"""
    
    @deprecated('filename')
    class c1(c):
        _loader_class = c.__mro__[2] #store the loader class for later use

        def __init__(self, filename:str=None, *args, **kwargs):
            if filename is not None:
                if not hasattr(self,'metadata'):
                    self.metadata = {}
                if hasattr(self,'_metadata_from_filename'):
                    self.metadata.update(self._metadata_from_filename(filename))
                self._loader_class.__init__(self, filename=os.path.abspath(filename), metadata=self.metadata)
            else:
                super().__init__(*args, **kwargs)

    #generate the docstring
    c1.__doc__ = c.__doc__
    c1._doc_params_ = {'Parameters':
                       """filename: str\n    Absolute or relative path to the file with model data. This argument is deprecated.""",
                       **c._doc_params_}
    c1.__init__.__doc__ = c1._generate_docstring()
    #update the call signature
    S = inspect.signature(c)
    S1 = inspect.signature(c1.__init__)
    #set default values to None if they are not set
    other_params = []
    for p in S.parameters.values():
        if p.default==p.empty:
            p = p.replace(default=None)
        other_params+=[p]
    params = [S1.parameters['self'],S1.parameters['filename'],*other_params]
    #fill the constructor signature
    c1.__init__.__signature__ = S.replace(parameters=params)
    #fill the class and module name to be the same as in class
    c1.__qualname__ = c.__qualname__
    c1.__name__ = c.__name__
    c1.__module__ = c.__module__
    #register the model in the list
    global all_models
    all_models.remove(c)
    all_models.add(c1)
    return c1
