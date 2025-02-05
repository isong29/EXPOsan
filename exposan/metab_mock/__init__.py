# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import os

folder = os.path.dirname(__file__)
data_path = os.path.join(folder, 'data')
results_path = os.path.join(folder, 'results')
figures_path = os.path.join(folder, 'figures')
# To save simulation results and generated figures
if not os.path.isdir(results_path): os.mkdir(results_path)
if not os.path.isdir(figures_path): os.mkdir(figures_path)

_impact_item_loaded = False
from exposan.metab import load_lca_data

from . import process
from .process import *

from . import units
from .units import *

from . import systems
from .systems import *


__all__ = (
    'folder',
    'data_path',
    'results_path',
    'figures_path',
    *process.__all__,
    *units.__all__,
    *systems.__all__,
	)


