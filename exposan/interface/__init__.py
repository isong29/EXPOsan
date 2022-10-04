# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os, qsdsan as qs

interface_path = os.path.dirname(__file__)
results_path = os.path.join(interface_path, 'results')
# To save simulation data
if not os.path.isdir(results_path): os.mkdir(results_path)


# %%

# =============================================================================
# Load the system
# =============================================================================

from . import system
from .system import *
_system_loaded = False
def _load_system():
    global sys, _system_loaded
    sys = create_system()
    _system_loaded = True


def load():
    if not _system_loaded: _load_system()
    dct = globals()
    dct.update(sys.flowsheet.to_dict())


def __getattr__(name):
    if not _system_loaded:
        raise AttributeError(f'module "{__name__}" not yet loaded, '
                             f'load module with `{__name__}.load()`.')


__all__ = (
    'interface_path',
    'results_path',
    *system.__all__,
)