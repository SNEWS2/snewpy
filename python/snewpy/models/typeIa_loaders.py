# -*- coding: utf-8 -*-
"""
The submodule ``snewpy.models.typeIa_loaders`` contains classes to load type Ia supernova
models from files stored on disk.
"""

# -*- coding: utf-8 -*-
"""
The submodule ``snewpy.models.pisn_loaders`` contains classes to load pair-instability supernova
models from files stored on disk.
"""

from snewpy.models.base import SNOwGLoBES

class Wright_2016(SNOwGLoBES):
    """The DDT model is described in 'Neutrinos from type Ia supernovae: The deflagration-to-detonation transition scenario', by Warren P. Wright et al.,
    [Phys. Rev. D94 (2016) 025026](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.94.025026), [arXiv:1605.01408](https://arxiv.org/abs/1605.01408).  
    There are 30 snapshots in time and the format of each data file is the SNOwGLoBES format.  

    The GCD model is described in 'Neutrinos from type Ia supernovae: The gravitationally confined detonation scenario' by Warren P. Wright et al.,
    [Phys. Rev. D95 (2017) 043006](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.95.043006), [arXiv:1609.07403](https://arxiv.org/abs/1609.07403).  
    There are 64 snapshots in time and the format of each data file is the SNOwGLoBES format.
    """
    pass





