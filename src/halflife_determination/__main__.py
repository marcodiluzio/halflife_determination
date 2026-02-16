# -*- coding: utf-8 -*-

"""
__main__.py of package halflife_determination
"""

import os
import sys
try:
    from halflife_determination import hl_elaboration as hle
except ImportError:
    import hl_elaboration as hle

print(hle.__doc__)

print(f"Arguments")
for i, arg in enumerate(sys.argv):
    if i > 0:
        print(f"Argument {i:>6}: {arg}")
    
try:
    path = sys.argv[1]
except IndexError:
    raise IndexError('argument defining the path of the folder containing the files is required')
if not os.path.isdir(path):
    raise FileNotFoundError('Folder not found!\ncheck the correct path or make sure path is a folder')

try:
    config_file = hle.load_config(sys.argv[2])
except (IndexError, TypeError, ValueError, FileNotFoundError, PermissionError):
    config_file = {}
    
defaults = {'apt':True, 'write_csv':True, 'MC_trials':10000, 'fit':'all', 'method':'all', 'output_path':'', 'iterative':False}
settings = {**defaults, **config_file}

half_life_results, additional_information = hle.elaboration(path, **settings)