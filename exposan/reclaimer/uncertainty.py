#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Hannah Lohman <hlohman94@gmail.com>
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

# Run uncertainty analysis and Spearman without country-specific settings
from exposan import reclaimer as re
from exposan.reclaimer import create_model, run_uncertainty

def run(model_IDs, seed=None, N=1000, country_specific=False, **model_kwargs):
    # Make it possible to run with one or more models
    if isinstance(model_IDs, str): model_IDs = (model_IDs, )
    for ID in model_IDs:
        model = create_model(ID, country_specific=country_specific, **model_kwargs)
        run_uncertainty(model, seed=seed, N=N)


if __name__ == '__main__':
    re.INCLUDE_RESOURCE_RECOVERY = True
    # run(('A', 'B', 'C', 'D'), seed=5, N=10)
    run(('B', 'C'), seed=5, N=50) # running systems B and C for contextual analysis with DMsan
