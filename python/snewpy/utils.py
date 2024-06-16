import numpy as np

def expand_dimensions_to(a:np.ndarray, ndim:int)->np.ndarray:
    """Expand the dimensions of the array, adding dimensions of len=1 to the right,
    so total dimensions equal to `ndim`"""
    new_shape = (list(a.shape)+[1]*ndim)[:ndim]
    return a.reshape(new_shape)