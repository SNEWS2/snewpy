from snewpy.flux import Flux
import numpy as np
import pytest

def test_init():
    data = np.random.random(size=(4,5,6))
    axes = {'x': range(4),
            'y': range(5),
            'z': range(6),
            }
    f = Flux(data, **axes)
    assert f.shape == data.shape
    assert np.allclose(f.array, data)
    assert f.axes.keys() == axes.keys()
    
def test_init_raises_on_inconsistent_axes():
    data = np.random.random(size=(4,5,6))
    with pytest.raises(ValueError):
        Flux(data, x=range(4), y=range(5), z=range(7))
    with pytest.raises(ValueError):
        Flux(data, x=range(4), y=range(5))
 
