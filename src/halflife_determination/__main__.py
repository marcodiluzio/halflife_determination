# -*- coding: utf-8 -*-

"""
__main__.py of package half-life_determination
"""

import os
import sys
from halflife_determination import hl_elaboration as hle

print(hle.__doc__)

print(f"Arguments")
for i, arg in enumerate(sys.argv):
    if i > 0:
        print(f"Argument {i:>6}: {arg}")
    
try:
    path = sys.argv[1]
except IndexError:
    path = os.path.abspath('.')

try:
    config_file = hle.load_config(sys.argv[2])
except (IndexError, TypeError, ValueError, FileNotFoundError, PermissionError):
    config_file = {}
    
defaults = {'apt':True, 'write_csv':True, 'MC_trials':10000, 'fit':'all', 'method':'all', 'output_path':'', 'iterative':False}
settings = {**defaults, **config_file}

half_life_results, additional_information = hle.elaboration(path, **settings)